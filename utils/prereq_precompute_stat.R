suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(arrow)
  library(jsonlite)
  library(stringr)
  library(tidyr)
})

# --- 1. Core Calculation Function (Your original code) ---
precompute_computeStat <- function(seurat_object){
  cat("Calculating summaries...\n")
  # Summary1
  Summary_perSubject_perStatus <- as.data.frame(seurat_object@meta.data %>%
                                                  group_by(.data$Subject,
                                                           .data$Status) %>%
                                                  tally())

  # Summary2 - by Clusters and Sample
  # Custom sort function
  custom_sort <- function(samples, suffix_order) {
    ids <- str_extract(samples, "^[^_]+_[^_]+")
    suffixes <- str_extract(samples, "[^_]+$")
    suffix_order <- match(suffixes, suffix_order)
    sorted <- samples[order(ids, suffix_order)]
  }
  suffix_order <- c("Prior", "Lesion", "Post")

  Summary_cluster_per_sample = seurat_object@meta.data %>%
    dplyr::count(.data$CellType_Level3, .data$orig.ident) %>%
    dplyr::rename(CellType = .data$CellType_Level3, 
                  Sample = .data$orig.ident) %>%
    tidyr::separate(.data$Sample, into = c("Subject", "Status"),
                    sep = "_", remove = FALSE) %>%
    mutate(
      Status = case_when(Status == "Entry" ~ "Prior",
                         Status == "8WPH" ~ "Post",
                         TRUE ~ Status)
    ) %>%
    mutate(
      CellType = factor(.data$CellType, levels = levels(seurat_object)),
      Sample = factor(.data$Sample,
                      levels = custom_sort(unique(.data$Sample), suffix_order)),
      Status = factor(.data$Status,
                      levels = c("Prior", "Lesion", "Post"))
    ) %>%
    dplyr::group_by(.data$Sample) %>%
    mutate(Total_Cells = sum(n)) %>%
    ungroup() %>%
    mutate(Percentage = n / .data$Total_Cells * 100) %>%
    dplyr::select(-.data$Total_Cells) %>%
    select(.data$Sample,
           .data$Subject,
           .data$Status,
           .data$CellType,
           n,
           .data$Percentage) %>%
    as.data.frame()

  # Summary3 by Cluster and Status
  n1 <- length(unique(seurat_object$Subject[seurat_object$Status == "Prior"]))
  n2 <- length(unique(seurat_object$Subject[seurat_object$Status == "Lesion"]))
  n3 <- length(unique(seurat_object$Subject[seurat_object$Status == "Post"]))

  Summary_cluster_per_status = seurat_object@meta.data %>%
    dplyr::count(.data$CellType_Level3, .data$Status) %>%
    dplyr::rename(CellType = .data$CellType_Level3) %>%
    mutate(
      CellType = factor(.data$CellType, levels = levels(seurat_object)),
      Status = factor(.data$Status, levels = c("Prior", "Lesion", "Post"))
    ) %>%
    dplyr::group_by(.data$Status) %>%
    mutate(Total_Cells = sum(n)) %>%
    ungroup() %>%
    mutate(Percentage = n / .data$Total_Cells * 100,
           Avg = n / case_when(
             Status == "Prior" ~ !!n1,
             Status == "Lesion" ~ !!n2,
             Status == "Post" ~ !!n3
           )) %>%
    dplyr::select(-.data$Total_Cells) %>%
    as.data.frame()

  cat("Summaries calculated.\n")
  return(list(Summary_perSubject_perStatus = Summary_perSubject_perStatus,
              Summary_cluster_per_sample = Summary_cluster_per_sample,
              Summary_cluster_per_status = Summary_cluster_per_status))
}

# --- 2. Wrapper Function for Loading, Computing, and Saving ---
process_and_save_stats <- function(seurat_rdata_path,
                                   seurat_object_name,
                                   output_prefix) {
  cat(paste("--- Processing stats for:", output_prefix, "---\n"))

  # Load the Seurat object
  cat(paste("Loading Seurat object from:", seurat_rdata_path, "\n"))
  load(seurat_rdata_path)
  seurat_obj <- get(seurat_object_name)
  cat(paste("Object", seurat_object_name, "loaded.\n"))

  # Compute the statistics
  stats_list <- precompute_computeStat(seurat_obj)

  # Save each summary to its own Parquet file
  cat("Saving results to Parquet files...\n")
  output_stats_subj_status_path <- file.path("DataWarehouse/Stat",
                                             paste0(output_prefix,
                                              "_stats_subject_status.parquet"))

  arrow::write_parquet(stats_list$Summary_perSubject_perStatus,
                       output_stats_subj_status_path)

  cat(paste("Saved Subject/Status stats to",
            output_stats_subj_status_path, "\n"))

  output_stats_clust_sample_path <- file.path("DataWarehouse/Stat",
                                              paste0(output_prefix,
                                                    "_stats_cluster_sample.parquet"))

  arrow::write_parquet(stats_list$Summary_cluster_per_sample,
                       output_stats_clust_sample_path)

  cat(paste("Saved Cluster/Sample stats to",
            output_stats_clust_sample_path, "\n"))

  output_stats_clust_status_path <- file.path("DataWarehouse/Stat",
                                              paste0(output_prefix,
                                              "_stats_cluster_status.parquet"))

  arrow::write_parquet(stats_list$Summary_cluster_per_status,
                       output_stats_clust_status_path)

  cat(paste("Saved Cluster/Status stats to",
            output_stats_clust_status_path, "\n"))

  cat(paste("--- Finished processing stats for", output_prefix, "---\n\n"))
}

# --- 3. Main Control Panel ---
# This block runs only when the script is executed from the command line
if (!interactive()) {
  cat("Starting pre-computation of statistical summaries...\n\n")

  # Settings for the Myeloid dataset
  process_and_save_stats(
    seurat_rdata_path = "DataLake/HSV434_integrated_Lev3_Myeloid.Rdata",
    seurat_object_name = "HSV434_integrated_Lev3_Myeloid",
    output_prefix = "myeloid"
  )

  # Settings for the T-cell dataset
  process_and_save_stats(
    seurat_rdata_path = "DataLake/HSV434_integrated_Lev3_Tcell.Rdata",
    seurat_object_name = "HSV434_integrated_Lev3_Tcell",
    output_prefix = "tcell"
  )

  cat("All statistical summaries have been pre-computed.\n")
}