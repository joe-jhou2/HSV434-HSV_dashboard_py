# utils/run_r_violin_plot.py

import base64
import os
import tempfile
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from utils.db_connection import dict_to_r_vector

def generate_violin_plot_from_df(df, color_file, selected_features):
    """
    Generate a violin plot using the R computeViolinPlot() function from a pandas DataFrame.
    """
    tmp_path = ""
    try:
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            tmp_path = tmp.name
        tmp_path_r = tmp_path.replace("\\", "/")

        r_color_string = dict_to_r_vector(color_file)

        # Pass the dataframe to R
        with localconverter(ro.default_converter + pandas2ri.converter):
            ro.globalenv["plot_data_filtered"] = df
            # ro.globalenv["color_file"] = color_file
            ro.globalenv['selected_features_r'] = ro.StrVector(selected_features)

        # color_file_r = color_file.replace("\\", "/")

        r_code = f"""
            suppressPackageStartupMessages({{
                library(ggplot2)
                library(dplyr)
                library(data.table)
                library(jsonlite) 
                library(svglite)
                library(ragg)
            }})
            
            # Suppressing warnings speeds up execution by preventing console I/O lag
            suppressWarnings(suppressMessages({{
            
                save_start <- Sys.time()
                cat("--- R Violin Plot: Script execution started ---\\n", file=stderr())
                
                CellType_color <- {r_color_string}
                if (is.list(CellType_color)) {{
                    CellType_color <- unlist(CellType_color)
                }}
                names(CellType_color) <- gsub(" ", "\n", names(CellType_color))

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

                # Bail out early if nothing to plot
                if (length(gene_cols) == 0) {{
                    p_out <- ggplot() + annotate("text", x=0.5, y=0.5, label="No data for filter selection") + theme_void()
                    ggsave("{tmp_path_r}", plot = p_out, width = 6, height = 4, dpi = 150)
                    stop("No data remaining after filtering.", call. = FALSE)
                }}

                # --- Generate Violin Plot ---
                exp_meta_df_reshape <- as.data.table(
                    plot_data_filtered[, !colnames(plot_data_filtered) %in% c("Subject", "Barcode", "UMAP_1", "UMAP_2"), drop = FALSE]
                )

                exp_meta_df_reshape <- data.table::melt(
                    exp_meta_df_reshape,
                    id.vars = c("CellType_Level3", "Status"),
                    variable.name = "Gene",
                    value.name = "Exp"
                )

                exp_meta_df_reshape$Exp = as.numeric(exp_meta_df_reshape$Exp)
                exp_meta_df_reshape$CellType_Level3 <- gsub(" ", "\n", exp_meta_df_reshape$CellType_Level3)
                
                # change gene name order (optional), prior data already done, but confirm
                exp_meta_df_reshape$CellType_Level3 = factor(exp_meta_df_reshape$CellType_Level3, levels = names(CellType_color))
                exp_meta_df_reshape$Status = factor(exp_meta_df_reshape$Status, levels = c("Prior", "Lesion", "Post"))
                exp_meta_df_reshape$Gene = factor(exp_meta_df_reshape$Gene, levels = gene_cols)
                exp_meta_df_reshape = exp_meta_df_reshape %>% droplevels()

                # plotting: each facet is gene by cluster, x is status, y is exp value,
                if (nrow(exp_meta_df_reshape) > 10000) {{
                exp_meta_df_reshape <- exp_meta_df_reshape[sample(.N, 10000)]
                message("Downsampled to 10k rows for faster violin rendering")
                }}

                g = ggplot(data = exp_meta_df_reshape,
                        aes(x = Status, y = Exp, fill = CellType_Level3)) +
                geom_violin(scale = 'width',
                            trim = TRUE,
                            draw_quantiles = c(0.25, 0.5, 0.75),
                            color = 'black',
                            size = 0.3,
                            alpha = 0.8) +
                stat_summary(fun = median, geom = "pointrange", color = "black", size = 0.2) +
                stat_summary(fun = median, geom = "line", color = "black", aes(group = 1)) +
                scale_fill_manual(values = CellType_color) +
                facet_grid(Gene ~ CellType_Level3, scales = "free_x") +
                theme_bw() +
                theme(
                    panel.spacing = unit(0.05, "lines"),
                    panel.grid = element_blank(),
                    strip.background = element_blank(),
                    strip.text.x = element_text(size = 15, angle = 0, face = "plain", hjust = 0.5),
                    strip.text.y = element_text(size = 15, angle = 0, face = "plain", hjust = 0.5),
                    axis.title.x = element_text(size = 0, color = 'black'),
                    axis.title.y = element_text(size = 15, color = 'black'),
                    axis.text.x = element_text(size = 15, color = 'black', angle = 90, hjust = 1, vjust = 0.5),
                    axis.text.y = element_text(size = 15, color = 'black'),
                    legend.position = "none"
                ) +
                labs(y = 'Log Normalized Expression')

                # --- Save plot safely ---
                n_genes = length(unique(exp_meta_df_reshape$Gene))
                n_cell_types = length(unique(exp_meta_df_reshape$CellType_Level3))

                width_factor = 2.0
                height_factor = 2.0
                    
                # Enforce minimums to prevent blur/stretch on small data
                final_w = max(4, n_cell_types * width_factor)
                final_h = max(3, n_genes * height_factor)

                ggsave("{tmp_path_r}", plot = g, device = "png", width = final_w, height = final_h, dpi = 300, limitsize = FALSE)

                end_time <- Sys.time()
                cat("Violin Plot Computation time: ", round(end_time - save_start, 3), " sec\n")
        }}) # End suppressWarnings
        )
        """
        ro.r(r_code)

        # Encode and return
        with open(tmp_path, "rb") as image_file:
            encoded_image = base64.b64encode(image_file.read()).decode()
        return f"data:image/png;base64,{encoded_image}", []

    except Exception as e:
        print(f"--- ERROR generating Violin Plot from df ---\\n{e}")
        return "/assets/error_placeholder.png", [str(e)]
    finally:
        if tmp_path and os.path.exists(tmp_path):
            try:
                os.remove(tmp_path)
            except Exception:
                pass
