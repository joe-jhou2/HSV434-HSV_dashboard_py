
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(arrow)
  library(jsonlite)
})

# Create a generic function to process a dataset ---
precompute_umap_extract <- function(seurat_rdata_path,
                                    seurat_object_name,
                                    color_rdata_path,
                                    color_object_name,
                                    output_prefix) {

  cat(paste("--- Precompute dataset for UMAP extraction:",
            output_prefix, "---\n"))

  # Load the Seurat and color objects
  load(seurat_rdata_path)
  load(color_rdata_path)

  # Get the objects by their string names
  seurat_obj <- get(seurat_object_name)
  color_obj <- get(color_object_name)

  # Extract UMAP coordinates and metadata using familiar R syntax
  cat("Extracting UMAP coordinates and metadata...\n")
  umap_coords <- as.data.frame(seurat_obj@reductions$umap@cell.embeddings) %>%
    rename_with(~"UMAP_1", matches("^umap_?1$", ignore.case = TRUE)) %>%
    rename_with(~"UMAP_2", matches("^umap_?2$", ignore.case = TRUE))

  metadata <- seurat_obj@meta.data

  # Combine into a single data frame for plotting
  plot_df <- cbind(metadata, umap_coords) %>%
    rownames_to_column("Barcode")

  # Define output paths
  output_data_path <- file.path("DataWarehouse/UMAP",
                                paste0(output_prefix, "_umap_data.parquet"))

  output_color_path <- file.path("DataWarehouse/Color",
                                 paste0(output_prefix, "_colors.json"))

  # Save the data using the arrow package (for Parquet)
  arrow::write_parquet(plot_df, output_data_path)
  cat(paste("✅ Saved UMAP data to", output_data_path, "\n"))

  # Save the color object using the jsonlite package
  color_map_as_list <- as.list(color_obj)
  jsonlite::write_json(color_map_as_list,
                       output_color_path,
                       auto_unbox = TRUE,
                       pretty = TRUE)

  cat(paste("✅ Saved color map to", output_color_path, "\n"))

  cat(paste("--- Finished processing", output_prefix, "---\n\n"))
}

# Define and run the pre-computation for all your datasets.
if (!interactive()) {

  # Settings for the Myeloid dataset
  precompute_umap_extract(
    seurat_rdata_path = "DataLake/HSV434_integrated_Lev3_Myeloid.Rdata",
    seurat_object_name = "HSV434_integrated_Lev3_Myeloid",
    color_rdata_path = "DataLake/HSV434_integrated_Lev3_Myeloid_color.Rdata",
    color_object_name = "HSV434_integrated_Lev3_Myeloid_color",
    output_prefix = "myeloid"
  )

  # Settings for the T-cell dataset
  precompute_umap_extract(
    seurat_rdata_path = "DataLake/HSV434_integrated_Lev3_Tcell.Rdata",
    seurat_object_name = "HSV434_integrated_Lev3_Tcell",
    color_rdata_path = "DataLake/HSV434_integrated_Lev3_Tcell_color.Rdata",
    color_object_name = "HSV434_integrated_Lev3_Tcell_color",
    output_prefix = "tcell"
  )

  cat("All datasets have been pre-computed.\n")
}