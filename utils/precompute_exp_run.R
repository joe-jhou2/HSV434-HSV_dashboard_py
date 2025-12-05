source("utils/precompute_exp.R")

if (!interactive()) {
  cat("Starting pre-computation of Expression data...\n\n")

  # --- Check if called dynamically (on-demand from Dash) ---
  extract_prefix <- Sys.getenv("EXTRACT_PREFIX", unset = NA)
  extract_genes  <- Sys.getenv("EXTRACT_GENES", unset = NA)

  if (!is.na(extract_prefix)) {
    # On-demand mode
    cat(paste("Dynamic mode: processing dataset", extract_prefix, "\n"))
    genes_to_fetch <- if (!is.na(extract_genes) && nzchar(extract_genes))
      strsplit(extract_genes, ",")[[1]] else NULL

    raw_data_base <- "DataLake"
    base_name <- paste0("HSV434_integrated_Lev3_",
                        tools::toTitleCase(extract_prefix))

    seurat_rdata <- file.path(raw_data_base, paste0(base_name, ".Rdata"))
    seurat_obj_name <- base_name

    process_exp_data(
      seurat_rdata_path = seurat_rdata,
      seurat_object_name = seurat_obj_name,
      output_prefix = extract_prefix,
      genes_to_extract = genes_to_fetch
    )

  } else {
    # --- Batch mode for full datasets ---
    core_marker_genes <- c("CD3D", "CD3E", "CD4", "CD8A",
                           "FOXP3", "IL2RA", "PDCD1", "CTLA4",
                           "IFNG", "GZMB", "TNF", "IL2")

    datasets_to_process <- list(
      list(prefix = "myeloid", genes = NULL),
      list(prefix = "tcell", genes = core_marker_genes)
    )

    raw_data_base <- "DataLake"
    obj_suffix <- ""

    for (ds in datasets_to_process) {
      prefix <- ds$prefix
      base_name <- paste0("HSV434_integrated_Lev3_", tools::toTitleCase(prefix))
      seurat_rdata <- file.path(raw_data_base, paste0(base_name, ".Rdata"))
      seurat_obj_name <- paste0(base_name, obj_suffix)

      if (!file.exists(seurat_rdata)) {
        cat(paste("❌ ERROR: Input RData not found for",
                  prefix, "at", seurat_rdata, "\n"))
        next
      }

      tryCatch({
        process_exp_data(
          seurat_rdata_path = seurat_rdata,
          seurat_object_name = seurat_obj_name,
          output_prefix = prefix,
          genes_to_extract = ds$genes
        )
      }, error = function(e) {
        cat(paste("❌ ERROR processing", prefix, ":", e$message, "\n"))
      })
    }
    cat("✅ All defined datasets processed.\n")
  }
}
