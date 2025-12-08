#!/usr/bin/env Rscript

# ============================================================
# Export all available gene names from a Seurat object to JSON
# Example usage:
# '''
# Rscript utils/export_available_genes.R DataLake/HSV434_integrated_Lev3_Tcell.Rdata  DataWarehouse/GEX/tcell_avail_genelist.json
# '''
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript utils/export_available_genes.R 
  <input_Rdata> <output_json>")
}

input_rdata <- args[1]
output_json <- args[2]

cat("Loading RData file:", input_rdata, "\n")
if (!file.exists(input_rdata)) {
  stop(paste("Input file not found:", input_rdata))
}

# --- Load Seurat object dynamically from RData ---
e <- new.env()
load(input_rdata, envir = e)

# The RData may contain multiple objects — find the Seurat one
seurat_obj_name <- names(e)[sapply(e, function(x) inherits(x, "Seurat"))]
if (length(seurat_obj_name) == 0) {
  stop("No Seurat object found in the loaded RData file.")
}

if (length(seurat_obj_name) > 1) {
  cat("Multiple Seurat objects found, using the first one:",
      seurat_obj_name[1], "\n")
}

seurat_obj <- e[[seurat_obj_name[1]]]

# --- Extract gene names ---
cat("Extracting gene names from RNA assay...\n")

if (!"RNA" %in% names(seurat_obj@assays)) {
  stop("RNA assay not found in the Seurat object.")
}

available_genes <- rownames(seurat_obj[["RNA"]])

cat(paste0("Found ", length(available_genes), " genes.\n"))

# --- Save to JSON ---
dir.create(dirname(output_json), recursive = TRUE, showWarnings = FALSE)
cat("Saving gene list to:", output_json, "\n")
write_json(as.list(available_genes),
           output_json, pretty = TRUE, auto_unbox = TRUE)

cat("Gene list successfully exported.\n")
