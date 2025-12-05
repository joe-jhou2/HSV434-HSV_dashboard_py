# --- 1. Load Required Libraries ---
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(arrow)
  library(jsonlite)
  library(stringr) # For matches()
  library(Matrix) # For handling sparse matrices, if needed
})

# Help function to calculate feature gene percentage
calculate_feature_pert <- function(feature_df,
                                   selectedfeature) {
  # Check for necessary columns in feature_df
  required_cols <- c("Subject", "CellType_Level3", "Status")

  if (!all(required_cols %in% names(feature_df))) {
    stop("Missing one or more: Subject, CellType_Level3, Status.")
  }

  # Ensure selectedfeature exist in feature_df
  if (!all(selectedfeature %in% names(feature_df))) {
    missing_genes <- selectedfeature[!selectedfeature %in% names(feature_df)]
    stop("Missing following genes: ",
         paste(missing_genes, collapse = ", "), ".")
  }

  # Process data
  # Step 1: Calculate positive percentages and total counts
  #         for each status and feature
  result_part1 <- feature_df %>%
    group_by(.data$Subject, .data$CellType_Level3, .data$Status) %>%
    summarise(
      Total_Cells = n(),  # Total cell count in each group
      across(all_of(selectedfeature), ~ sum(. > 0), .names = "n_{.col}"),
      .groups = "drop"
    )

  # Step 2: Calculate total positive cells across all statuses for each feature
  total_feature_plus <- result_part1 %>%
    group_by(.data$Subject, .data$CellType_Level3) %>%
    summarise(
      across(starts_with("n_"), sum, .names = "total_{.col}"),
      .groups = "drop"
    )

  # Step 3: Calculate the percentage contribution
  #.        for each status and final calculation
  final_result <- result_part1 %>%
    left_join(total_feature_plus, by = c("Subject", "CellType_Level3")) %>%
    mutate(across(
      starts_with("n_"),
      ~ . / get(paste0("total_", cur_column())) * . / Total_Cells * 100,
      .names = "{gsub('n_', '', .col)}"
    ))

  final_result <- final_result %>%
    select(.data$Subject, .data$CellType_Level3,
           .data$Status, all_of(selectedfeature))
}

# --- 2. Function to Process One Dataset ---
process_exp_data <- function(seurat_rdata_path,
                             seurat_object_name,
                             output_prefix,
                             genes_to_extract = NULL) {

  cat(paste("--- Processing expression data for:", output_prefix, "---\n"))

  # --- Load Data ---
  cat("   Loading RData files...\n")
  env <- new.env()
  load(seurat_rdata_path, envir = env)
  seurat_obj <- env[[seurat_object_name]]

  # Basic validation
  if (is.null(seurat_obj)) {
    stop("Seurat object '", seurat_object_name, "' not found.")
  }

  # --- Define Output Directory ---
  gex_dir <- file.path("DataWarehouse", "GEX")
  dir.create(gex_dir, recursive = TRUE, showWarnings = FALSE)

  pert_dir <- file.path("DataWarehouse", "Pert")
  dir.create(pert_dir, recursive = TRUE, showWarnings = FALSE)

  # --- Extract UMAP & Metadata ---
  cat("   Extracting UMAP coordinates and metadata...\n")
  umap_coords <- as.data.frame(seurat_obj@reductions$umap@cell.embeddings) %>%
    rename_with(~"UMAP_1", matches("^umap_?1$", ignore.case = TRUE)) %>%
    rename_with(~"UMAP_2", matches("^umap_?2$", ignore.case = TRUE))
  if (!all(c("UMAP_1", "UMAP_2") %in% colnames(umap_coords))) {
    stop("UMAP columns not found.")
  }

  metadata <- seurat_obj@meta.data %>% rownames_to_column("Barcode")

  # --- Extract RNA Expression Data ---
  expression_data_df <- NULL

  genes_to_fetch <- if (is.null(genes_to_extract)) {
    rownames(seurat_obj[["RNA"]])
  } else {
    genes_to_extract
  }

  action_desc <- if (is.null(genes_to_fetch)) {
    "ALL genes"
  } else {
    paste(length(genes_to_fetch), "specified genes")
  }

  cat(paste(
    "   Attempting to extract RNA expression data (",
    action_desc,
    ", data layer)...\n"
  ))

  tryCatch({
    if (!("RNA" %in% Assays(seurat_obj))) stop("RNA assay not found.")

    if (!is.null(genes_to_fetch)) {
      available_genes <- rownames(seurat_obj[["RNA"]])
      valid_genes <- intersect(genes_to_fetch, available_genes)
      invalid_genes <- setdiff(genes_to_fetch, available_genes)
      if (length(valid_genes) == 0) {
        stop(paste("None of requested genes found:",
                   paste(genes_to_fetch, collapse = ", ")))
      }

      if (length(invalid_genes) > 0) {
        cat(paste("   ⚠️ WARNING: Genes not found:",
                  paste(invalid_genes, collapse = ", "), "\n"))
      }
      genes_to_fetch <- valid_genes
    }

    # Extract data (features=NULL extracts all)
    rna_data <- FetchData(object = seurat_obj,
                          vars = genes_to_fetch,
                          layer = "data") %>%
      rownames_to_column("Barcode")

    expression_data_df <- rna_data

    cat(paste("   Successfully prepared expression data for",
              ncol(rna_data) - 1, "genes.\n"))

  }, error = function(e) {
    cat(paste("   ⚠️ ERROR during expression extraction:", e$message, "\n"))
    # If expression extraction fails, we proceed without it for this file
    cat("   Proceeding without expression data for this combined file.\n")
  })

  # --- Combine Metadata, UMAP, and (if available) Expression Data ---
  cat("   Combining extracted data...\n")

  combined_df <- metadata %>%
    left_join(umap_coords %>% rownames_to_column("Barcode"), by = "Barcode")

  # Join expression data if it was successfully extracted
  if (!is.null(expression_data_df)) {
    combined_df <- combined_df %>%
      left_join(expression_data_df, by = "Barcode")
  }

  gene_pert_df <- calculate_feature_pert(
    feature_df = combined_df,
    selectedfeature = genes_to_fetch
  )

  # --- Save to Parquet and Update Gene Index ---
  cat("   Saving data and updating index...\n")

  # Define output filename (timestamped if on-demand)
  extract_prefix <- Sys.getenv("EXTRACT_PREFIX", unset = "")
  is_dynamic_mode <- nzchar(extract_prefix)
  cat("   Dynamic mode:", is_dynamic_mode, "\n")

  if (is_dynamic_mode) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    gex_output_filename <- paste0(output_prefix,
                                  "_gex_",
                                  timestamp,
                                  ".parquet")
    pert_output_filename <- paste0(output_prefix,
                                   "_pert_",
                                   timestamp,
                                   ".parquet")
  } else {
    gex_output_filename <- paste0(output_prefix, "_gex_core.parquet")
    pert_output_filename <- paste0(output_prefix, "_pert_core.parquet")
  }
  gex_output_path <- file.path(gex_dir, gex_output_filename)
  pert_output_path <- file.path(pert_dir, pert_output_filename)

  # --- Save Combined Data File ---
  cat(paste("   Saving combined data to Parquet:", gex_output_path, "\n"))
  arrow::write_parquet(combined_df, gex_output_path)
  cat("   ✅ Saved combined data successfully.\n")

  # --- Save Gene Percentage Data File ---
  cat(paste("   Saving Gene Percentage data to Parquet:",
            pert_output_path, "\n"))
  arrow::write_parquet(gene_pert_df, pert_output_path)
  cat("   ✅ Saved Gene Percentage data successfully.\n")

  # --- Update Gene Index JSON ---
  if (!is.null(expression_data_df)) {
    cat("   Updating gene list index...\n")

    # Path to JSON
    output_gene_index_path <- file.path(gex_dir,
                                        paste0(output_prefix,
                                               "_gex_genes.json"))

    new_genes <- setdiff(colnames(expression_data_df), "Barcode")

    # Merge with existing index if available
    if (file.exists(output_gene_index_path)) {
      existing_genes <- jsonlite::read_json(output_gene_index_path,
                                            simplifyVector = TRUE)
      merged_genes <- sort(unique(c(existing_genes, new_genes)))
    } else {
      merged_genes <- sort(unique(new_genes))
    }

    jsonlite::write_json(merged_genes,
                         output_gene_index_path,
                         auto_unbox = TRUE)
    cat(paste("   ✅ Updated gene index to:", output_gene_index_path, "\n"))
  } else {
    cat("   ⚠️ No expression data extracted; index not updated.\n")
  }
  cat(paste("--- Finished processing:", output_prefix, "---\n\n"))
}