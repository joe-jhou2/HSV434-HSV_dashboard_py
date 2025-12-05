# utils/run_r_spatial_gene.py

import os
import base64
import re
import time
import boto3
from urllib.parse import urlparse
from dotenv import load_dotenv
from dash import html

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

# Load environment variables
load_dotenv()

# --- 1. DATA LOADING INFRASTRUCTURE ---
DATA_RELATIVE_PATH = "DataWarehouse/Visium/HSV434-Visium-A1_sub_label.Rdata"

def ensure_visium_data_exists():
    """
    Checks if the RData file exists locally. 
    If not, downloads it from S3 (ECS/New Environment support).
    Returns the R-compatible absolute path.
    """
    # Define paths
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    full_local_path = os.path.join(project_root, DATA_RELATIVE_PATH)

    # Check existence
    if not os.path.exists(full_local_path):
        print(f"📦 Visium Data missing: {DATA_RELATIVE_PATH}")
        print("🔄 Downloading from S3...")
        start_ts = time.time()

        try:
            # Create dir if missing
            os.makedirs(os.path.dirname(full_local_path), exist_ok=True)

            # Get Bucket Info
            bucket_uri = os.getenv("S3_BUCKET_URI")
            if not bucket_uri:
                raise ValueError("S3_BUCKET_URI is not set in .env")
            
            bucket_name = urlparse(bucket_uri).netloc
            
            # Construct S3 Key (Assumes repo structure matches S3 structure)
            s3_key = f"Joe/HSV_Dashboard_py/{DATA_RELATIVE_PATH}"

            # Download
            s3 = boto3.client('s3')
            s3.download_file(bucket_name, s3_key, full_local_path)
            
            print(f"✅ Download complete in {time.time() - start_ts:.2f}s")

        except Exception as e:
            # If we can't get the data, the app can't function. Raise critical error.
            raise FileNotFoundError(f"CRITICAL: Could not download Visium RData from S3. {e}")

    # Return path with forward slashes (Required for R on Windows/Linux)
    return full_local_path.replace("\\", "/")

# --- 2. INITIALIZE R ENVIRONMENT ---

# Get the valid data path (Trigger download if needed)
try:
    R_DATA_PATH = ensure_visium_data_exists()
except Exception as e:
    print(f"❌ Error initializing Spatial Data: {e}")
    R_DATA_PATH = None

# Load necessary R libraries
ro.r('''
suppressMessages({
  suppressPackageStartupMessages({
    suppressWarnings({
      library(Seurat)
      library(SeuratObject)
      library(sp)
      library(ggplot2)
      library(cowplot)
    })
  })
})
''')

# Load the Visium dataset
if R_DATA_PATH:
    print(f"R is loading data from: {R_DATA_PATH}")

    ro.r('''
      suppressMessages({
        suppressWarnings({
          load("DataWarehouse/Visium/HSV434-Visium-A1_sub_label.Rdata")
        })
      })
    ''')

    # Get gene names from the SCT assay's data slot
    original_gene_names = ro.r('rownames(A1_sub@assays$SCT$data)')

    # Create the mapping from an uppercase gene name to its original, correct case
    UPPER_TO_ORIGINAL_CASE_MAP = {name.upper(): name for name in original_gene_names}

    # Create the set of uppercase names for fast, case-insensitive validation
    AVAILABLE_GENES_UPPER = set(UPPER_TO_ORIGINAL_CASE_MAP.keys())
    # print(f"Validation enabled. {len(AVAILABLE_GENES_UPPER)} genes available from SCT assay.")

    # Define the theme once in the global R environment
    ro.r("""
        spatial_feature_theme <- theme(
            plot.background = element_rect(fill = "white"),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            panel.background = element_rect(fill = "white"),
            legend.position = "right",
            legend.margin = margin(t = 0.1),
            plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 8)
        )
    """)
else:
    # Fallback if data failed to load (prevents import crash, but features will fail)
    AVAILABLE_GENES_UPPER = set()
    UPPER_TO_ORIGINAL_CASE_MAP = {}

# -----------------------------------------------------------------------------
def run_r_spatial_gene(gene_names_str):
    """
    Validates gene names, then runs R SpatialFeaturePlot for valid genes.
    Handles single or multiple comma-separated genes.
    """
    if not gene_names_str:
        return html.Div("⚠️ Please enter at least one gene name.", style={"color": "red"}), False

    # Make input case-insensitive and validate against available genes
    # input_genes = [gene for gene in re.split(r'[\s,;]+', gene_names_str) if gene]
    # input_genes_upper = {gene.strip().upper() for gene in input_genes}
    # valid_genes_upper = sorted(list(input_genes_upper.intersection(AVAILABLE_GENES_UPPER)))
    # invalid_genes_upper = sorted(list(input_genes_upper.difference(AVAILABLE_GENES_UPPER)))
    # Split & clean
    input_genes = [gene.strip() for gene in re.split(r'[\s,;]+', gene_names_str) if gene.strip()]

    # Preserve input order, validate case-insensitively
    valid_genes_original_case = []
    invalid_genes_upper = []

    for gene in input_genes:
        g_up = gene.upper()
        if g_up in AVAILABLE_GENES_UPPER:
            valid_genes_original_case.append(UPPER_TO_ORIGINAL_CASE_MAP[g_up])
        else:
            invalid_genes_upper.append(g_up)

    if not valid_genes_original_case:
        error_msg = f"❌ Error: None of the entered genes were found. Invalid: {', '.join(invalid_genes_upper)}"
        return html.Pre(error_msg), False

    # Convert valid genes back to their ORIGINAL case for R
    r_gene_vector = 'c({})'.format(','.join([f'"{gene}"' for gene in valid_genes_original_case]))

    num_rows = -(-len(valid_genes_original_case) // 4)
    plot_height = max(4, 4 * num_rows)
    plot_width = 16

    try:
        r_code = f"""

        genes_to_check <- {r_gene_vector}
        
        # Check for non-zero expression in the SCT data slot
        sct_data <- GetAssayData(A1_sub, assay = "SCT", layer = "data")[genes_to_check, , drop = FALSE]
        plottable_genes <- names(which(rowSums(as.matrix(sct_data)) > 0))
        
        if (length(plottable_genes) == 0) {{
          "NO_PLOTTABLE_GENES"
        }} else {{
          plot_list <- lapply(plottable_genes, function(gene) {{
              p <- SpatialFeaturePlot(A1_sub, 
                                      features = gene, 
                                      pt.size.factor = 4,
                                      image.alpha = 0, 
                                      stroke = NA, 
                                      alpha = 1, 
                                      max.cutoff = 'q95') +
                  spatial_feature_theme
              return(p)
          }})
          combined_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = 4)
          tmp_path <- tempfile(fileext=".png")
          ggsave(tmp_path, plot = combined_plot, width = {plot_width}, height = {plot_height}, units = "in", dpi = 150)
          tmp_path
        }}
        """

        with localconverter(ro.default_converter + pandas2ri.converter):
            output_path = ro.r(r_code)[0]

        # Check for the special message from R
        if output_path == "NO_PLOTTABLE_GENES":
            unplottable = ', '.join(set(valid_genes_original_case))
            error_msg = f"❌ Error: The gene(s) '{unplottable}' exist but have no expression data to plot."
            return html.Pre(error_msg), False

        # proceed with creating the image
        with open(output_path, "rb") as f:
            encoded = base64.b64encode(f.read()).decode()
        os.remove(output_path)
        
        plot_img = html.Img(src=f"data:image/png;base64,{encoded}",
                            style={"width": "100%", "border": "1px solid #ccc"})

        # Show a warning for invalid ones
        if invalid_genes_upper:
            warning_msg = html.Div([
                html.Strong("Warning:"),
                f" The following genes were not found and were ignored: {', '.join(invalid_genes_upper)}"
            ], style={'padding': '10px', 'backgroundColor': '#FFF3CD', 'color': '#664D03', 'border': '1px solid #FFECB5', 'borderRadius': '5px', 'marginBottom': '15px'})
            return html.Div([warning_msg, plot_img]), False
        else:
            return plot_img, False

    except Exception:
        import traceback
        error_details = traceback.format_exc()
        return html.Pre(f"❌ An error occurred during R execution:\n\n{error_details}"), False