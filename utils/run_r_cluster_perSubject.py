import base64
import os
import tempfile
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

def generate_PerSubject_StackBar_plots(dataset_prefix, subjects=None):
    """
    Generates a stacked bar plot showing cell type proportions per subject, faceted by status.
    Accepts an optional list of subjects for filtering. Returns a Base64 image string.
    """
    # Define necessary file paths
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    stats_file = os.path.join(project_root, "DataWarehouse/Stat", f"{dataset_prefix}_stats_cluster_sample.parquet")
    color_file = os.path.join(project_root, "DataWarehouse/Color", f"{dataset_prefix}_colors.json")

    # --- Convert selected subjects to an R vector string ---
    subjects_r_vector = "NULL"
    if subjects:
        subjects_r_vector = f"c({', '.join([f'{chr(39)}{s}{chr(39)}' for s in subjects])})"

    # --- Calculate dynamic height based on number of subjects ---
    # Create a secure, temporary file for the R plot
    tmp_path = ""
    try:
        # Check if required files exist
        if not os.path.exists(stats_file):
            raise FileNotFoundError(f"Stats file not found: {stats_file}")
        if not os.path.exists(color_file):
            raise FileNotFoundError(f"Color file not found: {color_file}")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            tmp_path = tmp.name

        # Construct the self-contained R script as a string
        r_code_string = f"""
            suppressPackageStartupMessages({{
              library(ggplot2)
              library(dplyr)
              library(arrow)
              library(jsonlite)
              library(tidyr)    # Needed for complete()
              library(forcats)  # Needed for fct_rev()
            }})

            # Read the pre-computed statistics and colors
            Summary_cluster_per_sample <- arrow::read_parquet("{stats_file}")
            CellType_color <- jsonlite::fromJSON("{color_file}")
            selectedSubject <- {subjects_r_vector} # Assign the vector from Python

            # Define known subject levels for consistent ordering (adjust if needed)
            subject_levels <- c("Subject1", "Subject2", "Subject3", "Subject5", "Subject6",
                                "Subject7", "Subject8", "Subject9", "Subject10", "Subject11","Subject12",
                                "Subject13", "Subject14", "Subject15", "Subject16", "Subject17", "Subject18")
            status_levels <- c("Prior", "Lesion", "Post")

            # --- Data Wrangling (adapted from your R code) ---
            expanded_data <- Summary_cluster_per_sample %>%
              # Ensure Percentage column exists and is numeric, n might be needed too
              mutate(Percentage = as.numeric(Percentage), n = as.numeric(n)) %>%
              # Complete combinations to ensure all bars appear
              tidyr::complete(CellType, Subject, Status, fill = list(Percentage = 0, n = 0)) %>%
              # Apply subject filtering
              dplyr::filter(if (is.null(selectedSubject)) TRUE else Subject %in% selectedSubject) %>%
              # Set factor levels for plotting order
              mutate(
                  Subject = fct_rev(factor(Subject, levels = intersect(subject_levels, unique(.$Subject)))), # Use only relevant levels
                  Status = factor(Status, levels = status_levels),
                  CellType = factor(CellType) # Let ggplot handle CellType order or define levels if needed
              )

            # Re-estimate height if 'All' was selected initially
            actual_num_subjects <- n_distinct(expanded_data$Subject)
            plot_height_final <- 1.3 + actual_num_subjects * 0.3

            # Check if data exists after filtering
            if (nrow(expanded_data) == 0) {{
                p <- ggplot() + annotate("text", x=0.5, y=0.5, label="No data for selection") + theme_void()
            }} else {{
                p <- ggplot(expanded_data, aes(x = Subject, y = Percentage / 100, fill = CellType)) + # Divide Percentage by 100 for scale_y_continuous
                  geom_col(position = position_fill(reverse = TRUE), color = "black", linewidth = 0.2) + # Use linewidth
                  scale_y_continuous(labels = scales::percent) +
                  xlab("") +
                  ylab("Proportion") +
                  facet_grid(. ~ Status) +
                  scale_fill_manual(values = CellType_color) +
                  theme_bw(base_size = 10) + # Adjust base size for readability
                  theme(
                      panel.border = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.spacing.x = unit(1, "lines"), # Adjust horizontal spacing between facets
                      strip.background = element_rect(fill="grey90", color=NA), # Style facet titles
                      strip.text = element_text(size = 12, face = "bold"),
                      axis.text = element_text(size = 9, color = "black"), # Adjust text size
                      axis.title = element_text(size = 11, color = "black"),
                      axis.ticks.y = element_blank(), # Remove y-axis ticks with coord_flip
                      axis.line.x = element_line(color = "black"),
                      legend.position = "bottom",
                      legend.title = element_blank(),
                      legend.text = element_text(size = 9),
                      legend.key.size = unit(0.4, "cm") # Adjust legend key size
                  ) +
                  guides(fill = guide_legend(ncol = 8, byrow = TRUE)) +
                  coord_flip() # Flip coordinates
            }}

            # Save the plot using calculated height
            ggsave("{tmp_path}", plot = p, width = 10, height = plot_height_final, dpi = 200)
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
        print(f"--- ERROR generating Subject StackBar plot for dataset: {dataset_prefix} ---")
        print(traceback.format_exc())
        return "/assets/error_placeholder.png" # Return placeholder on error
    finally:
        # Ensure the temporary file is always removed
        if os.path.exists(tmp_path):
            os.remove(tmp_path)