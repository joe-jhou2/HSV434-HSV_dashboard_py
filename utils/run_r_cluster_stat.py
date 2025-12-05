import base64
import os
import tempfile
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

def generate_clusterStat_plots(dataset_prefix):
    """
    Generates combined count and percentage bar plots for cell types vs. status
    by executing a self-contained R script. Returns a Base64 image string.
    """
    # Define necessary file paths
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    # Only need the cluster_status summary file for these plots
    stats_file = os.path.join(project_root, "DataWarehouse/Stat", f"{dataset_prefix}_stats_cluster_status.parquet")

    # Create a secure, temporary file for the R plot
    tmp_path = ""
    try:
        # Check if the required stats file exists first
        if not os.path.exists(stats_file):
            raise FileNotFoundError(f"Statistics file not found: {stats_file}. Please run pre-computation.")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            tmp_path = tmp.name

        # Construct the self-contained R script as a string
        r_code_string = f"""
            suppressPackageStartupMessages({{
              library(ggplot2)
              library(dplyr)
              library(arrow)
              library(patchwork) # Load patchwork for combining plots
            }})

            # Read the pre-computed statistics data
            Summary_cluster_per_status <- arrow::read_parquet("{stats_file}")

            # Define the shared theme
            theme_BarCellType2Status = theme(
                panel.background = element_rect(color = "black", fill = NA),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                strip.text = element_text(size=10, face="bold"), # Show strip text for facets
                strip.background = element_rect(fill="grey90", color="black"),
                plot.title = element_text(size = 14, hjust = 0.5, face="bold"),
                axis.title.x = element_blank(),
                axis.title.y = element_text(size = 12, color = "black"),
                axis.text.x = element_text(size = 10, color = "black", angle = 90, hjust = 1, vjust = 0.5), # Angle text for better fit
                axis.text.y = element_text(size = 10, color = "black"),
                legend.title = element_blank(),
                legend.text = element_text(size = 11, color = "black"),
                legend.position = "bottom"
            )

            # Count plot (p1)
            p1 = ggplot(Summary_cluster_per_status, aes(x = Status, y = n, fill = Status)) +
              geom_bar(stat = "identity", position = position_dodge(width = 1)) +
              scale_fill_manual(values = c("Prior" = "#0F9D58", "Lesion" = "#DB4437", "Post" =  "#F4B400")) +
              ylab("Count") +
              ggtitle("Counts per Status") +
              facet_wrap(~ CellType, scales = "free_y", ncol = 6) + # Facet by CellType, free y-axis
              theme_bw(base_size = 10) + # Use theme_bw and adjust base size
              theme_BarCellType2Status

            # Percentage plot (p2)
            p2 = ggplot(Summary_cluster_per_status, aes(x = Status, y = Percentage, fill = Status)) +
              geom_bar(stat = "identity", position = position_dodge(width = 1)) +
              scale_fill_manual(values = c("Prior" = "#0F9D58", "Lesion" = "#DB4437", "Post" =  "#F4B400")) +
              ylab("Percentage (%)") +
              ggtitle("Percentage per Status") +
              facet_wrap(~ CellType, scales = "free_y", ncol = 6) + # Facet by CellType, free y-axis
              theme_bw(base_size = 10) +
              theme_BarCellType2Status

            # Combine the plots using patchwork (p1 on top of p2)
            combined_plot = p1 / p2

            # Save the combined plot
            n_cell_types = length(unique(Summary_cluster_per_status$CellType))
            plot_height = max(10, n_cell_types * 0.55)
            ggsave("{tmp_path}", plot = combined_plot, width = 10, height = plot_height, dpi = 200)
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
        print(f"--- ERROR generating cluster statistics plots for dataset: {dataset_prefix} ---")
        print(traceback.format_exc())
        # Return path to an error placeholder image if you have one
        return "/assets/error_placeholder.png"
    finally:
        # Ensure the temporary file is always removed
        if os.path.exists(tmp_path):
            os.remove(tmp_path)