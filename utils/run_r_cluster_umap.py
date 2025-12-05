# utils/run_r_cluster_umap.py

import base64
import os
import tempfile
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

def generate_umap_plot(dataset_prefix, status="All", title="", clusters=None, subjects=None):
    """
    Generates a ggplot UMAP by executing a self-contained R script that
    saves to a temp file, which Python then reads and encodes.
    """
    # Define all necessary file paths
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    data_file = os.path.join(project_root, "DataWarehouse/UMAP", f"{dataset_prefix}_umap_data.parquet")
    color_file = os.path.join(project_root, "DataWarehouse/Color", f"{dataset_prefix}_colors.json")

    clusters_r_vector = "NULL"
    if clusters:
        clusters_r_vector = f"c({', '.join([f'{chr(39)}{c}{chr(39)}' for c in clusters])})"

    subjects_r_vector = "NULL"
    if subjects:
        subjects_r_vector = f"c({', '.join([f'{chr(39)}{s}{chr(39)}' for s in subjects])})"

    legend_r_code = ''
    if status != "All":
        legend_r_code = 'p <- p + theme(legend.position = "none")'

    #  Create a secure, temporary file for the R plot
    tmp_path = ""
    try:
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            tmp_path = tmp.name

        # Construct the self-contained R script as a string
        r_code_string = f"""
            suppressPackageStartupMessages({{
              library(ggplot2)
              library(dplyr)
              library(arrow)
              library(jsonlite)
            }})

            # Read the data and color files
            plot_df <- arrow::read_parquet("{data_file}")
            cell_colors <- jsonlite::fromJSON("{color_file}")

            # Assign the filter vectors from Python
            selected_clusters <- {clusters_r_vector}
            selected_subjects <- {subjects_r_vector}

            # Filter data based on status
            if ("{status}" != "All") {{
                plot_df <- plot_df %>% filter(Status == "{status}")
            }}

            # Filter data based on selected clusters
            if (!is.null(selected_clusters)) {{
                plot_df <- plot_df %>% filter(CellType_Level3 %in% selected_clusters)
            }}

            # Filter data based on selected subjects
            if (!is.null(selected_subjects)) {{
                plot_df <- plot_df %>% filter(Subject %in% selected_subjects)
            }}


            if (nrow(plot_df) == 0) {{
                p <- ggplot() + annotate("text", x=0, y=0, label="No data for selection") + theme_void()
            }} else {{
                scale_limit <- ceiling(max(abs(plot_df[, c("UMAP_1", "UMAP_2")])))
                
                p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
                    geom_point(aes(color = CellType_Level3), size = 0.5, alpha = 0.8) +
                    scale_color_manual(values = cell_colors) +
                    labs(x = "UMAP1", y = "UMAP2", title = "{title}", color = NULL) +
                    coord_fixed(ratio = 1) + 
                    xlim(-scale_limit, scale_limit) + 
                    ylim(-scale_limit, scale_limit) +
                    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
                          panel.background = element_rect(fill = NA), 
                          panel.grid = element_blank(),
                          axis.title = element_text(size = 14), 
                          legend.title = NULL,
                          legend.text = element_text(size = 10),
                          legend.spacing.x = unit(0.01, 'cm'),
                          legend.spacing.y = unit(0.01, 'cm'),
                          legend.position = "right",
                          legend.key.size = unit(0.5, 'cm'),
                          legend.box.margin = margin(0, 0, 0, 0))+
                   guides(color = guide_legend(override.aes = list(size = 2), title.hjust = 0.5, ncol = 2, byrow = TRUE, direction = "horizontal"))
            }}
            
            {legend_r_code}

            # Save the plot to the temporary file path provided by Python
            ggsave("{tmp_path}", plot = p, width = 7, height = 5, dpi = 200)
        """

        # Execute the R code
        with localconverter(ro.default_converter + pandas2ri.converter):
            ro.r(r_code_string)

        # Read the generated image file and encode it
        with open(tmp_path, "rb") as image_file:
            encoded_image = base64.b64encode(image_file.read()).decode()
        
        return f"data:image/png;base64,{encoded_image}"

    except Exception:
        import traceback
        print(f"--- ERROR generating UMAP for Status: {status} ---")
        print(traceback.format_exc())
        return ""
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)