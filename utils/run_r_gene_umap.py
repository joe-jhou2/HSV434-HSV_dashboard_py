# utils/run_r_gene_umap.py

import base64
import os
import tempfile
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

def generate_feature_umap_from_df(df, selected_features):
    """
    Generate UMAP plots (colored by gene expression or clusters)
    directly from a pre-loaded pandas DataFrame.
    """
    tmp_path = ""
    try:
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            tmp_path = tmp.name
        tmp_path_r = tmp_path.replace("\\", "/")

        # Pass the dataframe to R
        with localconverter(ro.default_converter + pandas2ri.converter):
            ro.globalenv["plot_data_filtered"] = df
            ro.globalenv['selected_features_r'] = ro.StrVector(selected_features)

        r_code = f"""
            suppressPackageStartupMessages({{
                library(ggplot2)
                library(dplyr)
                library(patchwork)
                library(scales)
                library(grid)
            }})

            # Suppressing warnings speeds up execution by preventing console I/O lag
            suppressWarnings(suppressMessages({{
            
                # Start timing
                save_start <- Sys.time()
                cat("--- R Feature UMAP Plot: Script execution started ---\\n", file=stderr())

                # --- Determine genes available for plotting ---
                gene_cols <- setdiff(colnames(plot_data_filtered),
                                    c("Barcode", "UMAP_1", "UMAP_2", "CellType_Level3", "Subject", "Status"))
                cat(paste0("   Genes available for plotting: ", paste(gene_cols, collapse=", "), "\\n"), file=stderr())

                selectedFeature = selected_features_r
                if (length(selectedFeature) > 0) {{
                    gene_cols <- intersect(selectedFeature, gene_cols)
                }}

                # reorder genes to match selectedFeature order
                gene_cols <- gene_cols[match(gene_cols, selectedFeature)]

                if (length(gene_cols) > 0) {{
                    # base UMAP
                    scale_limit <- ceiling(max(abs(plot_data_filtered[, c("UMAP_1", "UMAP_2")]), na.rm=TRUE))
                    base_umap_plot <- ggplot(data = plot_data_filtered, aes(x = UMAP_1, y = UMAP_2)) +
                        geom_density_2d(aes(color = ..level..), size = 0.5, bins = 15, color="grey70") +
                        scale_color_identity() +
                        coord_fixed(ratio = 1) +
                        xlim(-scale_limit, scale_limit) + ylim(-scale_limit, scale_limit) +
                        theme_bw() +
                        theme(panel.grid = element_blank(),
                            axis.title = element_text(size = 10),
                            legend.position = "none")

                    plot_list <- lapply(gene_cols, function(singleFeature) {{
                        expressing_cells <- plot_data_filtered %>%
                            filter(!is.na(.data[[singleFeature]]) & .data[[singleFeature]] > 0)

                        max_expr <- max(expressing_cells[[singleFeature]], 0, na.rm = TRUE)
                        if (is.na(max_expr) || max_expr <= 0) max_expr <- 1

                        gg <- base_umap_plot +
                            geom_point(data = expressing_cells,
                                    aes(color = .data[[singleFeature]]),
                                    size = 0.2, alpha = 1) +
                            scale_color_gradient(low = "lightyellow", high = "red",
                                                limits = c(0, max_expr), oob = scales::squish) +
                            labs(title = singleFeature, color = "Expr.") +
                            theme(legend.position = "right",
                                legend.title = element_text(size = 9),
                                legend.text = element_text(size = 8),
                                legend.key.size = unit(0.4, 'cm'),
                                plot.title = element_text(hjust = 0.5, face = "bold"))
                        return(gg)
                    }})
                    final_plot <- patchwork::wrap_plots(plot_list, ncol = min(4, length(plot_list)))
                }} else {{
                    # fallback: cluster coloring
                    scale_limit <- ceiling(max(abs(plot_data_filtered[, c("UMAP_1", "UMAP_2")]), na.rm=TRUE))
                    final_plot <- ggplot(data = plot_data_filtered, aes(x = UMAP_1, y = UMAP_2)) +
                        geom_point(aes(color = CellType_Level3), size = 0.4, alpha = 0.8) +
                        labs(x = "UMAP1", y = "UMAP2", title = "Clusters", color = "CellType") +
                        coord_fixed(ratio = 1) +
                        xlim(-scale_limit, scale_limit) + ylim(-scale_limit, scale_limit) +
                        theme_bw() +
                        guides(colour = guide_legend(override.aes = list(size = 3)))
                }}

                # --- Save plot ---
                plot_width <- 10
                plot_height <- if (length(gene_cols) > 0) max(4, 2.5 * ceiling(length(gene_cols) / min(4, length(gene_cols)))) else 6

                ggsave("{tmp_path_r}", plot = final_plot,
                    width = plot_width, height = plot_height,
                    dpi = 150, limitsize = FALSE)
                
                # End timing
                end_time <- Sys.time()
                cat("Feature UMAP Computation time: ", round(end_time - save_start, 3), " sec\n")
            }}))
        """
        ro.r(r_code)

        # Encode and return
        with open(tmp_path, "rb") as image_file:
            encoded_image = base64.b64encode(image_file.read()).decode()
        return f"data:image/png;base64,{encoded_image}", []

    except Exception as e:
        print(f"--- ERROR generating feature UMAP from df ---\\n{e}")
        return "/assets/error_placeholder.png", [str(e)]
    finally:
        if tmp_path and os.path.exists(tmp_path):
            try:
                os.remove(tmp_path)
            except Exception:
                pass
