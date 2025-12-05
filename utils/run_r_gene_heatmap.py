# utils/run_r_heatmap.py

import base64
import os
import tempfile
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from utils.db_connection import dict_to_r_vector

def generate_heatmap_from_df(df, color_file, selected_features):
    """
    Generate a heatmap using the R computeHeatmap() function from a pre-loaded DataFrame.
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
                library(arrow)
                library(duckdb)
                library(glue)
                library(jsonlite)
                library(RColorBrewer)
                library(circlize)
                library(ComplexHeatmap)
                library(ragg)
            }})

            # Suppressing warnings speeds up execution by preventing console I/O lag
            suppressWarnings({{
                # Start timing
                start_time <- Sys.time()
                cat("--- R Heatmap: Script execution started ---\\n", file=stderr())
                
                CellType_color <- {r_color_string}

                if (is.list(CellType_color)) {{
                    CellType_color <- unlist(CellType_color)
                }}

                # Bail out early if nothing to plot
                if (nrow(plot_data_filtered) == 0) {{
                    p_out <- ggplot() + annotate("text", x=0.5, y=0.5, label="No data for filter selection") + theme_void()
                    ggsave("{tmp_path_r}", plot = p_out, width = 6, height = 4, dpi = 150)
                    stop("No data remaining after filtering.", call. = FALSE)
                }}

                selectedFeature = selected_features_r
                if (length(selectedFeature) == 0) {{
                    stop("❌ No selected features provided from Python.")
                }}

                # --- Generate Heatmap ---
                # Prepare matrix for heatmap
                ht_data = as.matrix(t(plot_data_filtered[, !(names(plot_data_filtered) %in% c("orig.ident", "Barcode", "UMAP_1", "UMAP_2", "Subject", "CellType_Level3", "Status"))]))

                # Prepare annotation for heatmap
                ht_meta = plot_data_filtered[, c("CellType_Level3", "Status", "Subject"), drop = FALSE]
                ht_meta$CellType_Level3 = factor(ht_meta$CellType_Level3, levels = names(CellType_color))
                ht_meta$Status = factor(ht_meta$Status, levels = c("Prior", "Lesion", "Post"))

                # re-order columns by celltype --> status
                column_order = order(plot_data_filtered$CellType_Level3, plot_data_filtered$Status, plot_data_filtered$Subject)
                
                # re-order rows by selected features
                row_order = match(selectedFeature, rownames(ht_data))
                ht_data = ht_data[row_order, ]
                
                # Create column annotation for CellType
                anno_celltype_status = HeatmapAnnotation(CellType = ht_meta$CellType_Level3,
                                                        Status = ht_meta$Status,
                                                        annotation_name_gp = gpar(fontsize = 12), 
                                                        border = TRUE,
                                                        simple_anno_size = unit(0.4, "cm"), 
                                                        annotation_legend_param = list(title_gp = gpar(fontsize = 12), 
                                                                                        labels_gp = gpar(fontsize = 12),
                                                                                        nrow = 3, by_row = FALSE, direction = "horizontal"),
                                                        col = list(CellType = CellType_color,
                                                                    Status = c("Prior" = "#0F9D58",  "Lesion" = "#DB4437", "Post" = "#F4B400")))
                
                # Create a color mapping for the heatmap
                # c(min(ht_data), max(ht_data))
                col_fun = colorRamp2(c(0, 0.01, 7), c("gray","blue", "red"))
                ht_opt$message = FALSE
                
                # Generate the heatmap with column annotations
                htmp = Heatmap(ht_data, col = col_fun, row_names_side = "right", 
                            cluster_rows = FALSE, cluster_columns = FALSE, 
                            show_column_names = FALSE, column_order = column_order, 
                            top_annotation = anno_celltype_status,
                            row_split = seq(1, nrow(ht_data)),
                            column_title_gp = gpar(fontsize = 0), 
                            column_split = ht_meta$CellType_Level3, 
                            column_title_rot = 45, 
                            row_names_gp = gpar(fontsize = 12),
                            row_title = NULL,
                            row_gap = unit(1, "mm"),
                            border = TRUE, 
                            heatmap_legend_param = list(title = "Expression\nLevel", title_position = "topcenter", direction = "horizontal", 
                                                            title_gp = gpar(col = "red", fontsize = 12), 
                                                            labels_gp = gpar(col = "black", fontsize = 12)), use_raster = FALSE)

                # --- Save plot safely ---
                save_start <- Sys.time()
                cat("  Step 4: Saving ggplot heatmap...\n", file=stderr())

                agg_png("{tmp_path_r}", width = 12, height = max(3, 0.4*dim(ht_data)[1]), units = "in", res = 150)
                draw(htmp, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
                dev.off()

                # End timing
                end_time <- Sys.time()
                cat("Heatmap Computation time: ", round(end_time - save_start, 3), " sec\n")
            }})
        """
        ro.r(r_code)

        # Encode and return
        with open(tmp_path, "rb") as image_file:
            encoded_image = base64.b64encode(image_file.read()).decode()
        return f"data:image/png;base64,{encoded_image}", []

    except Exception as e:
        print(f"--- ERROR generating Heatmap Plot from df ---\\n{e}")
        return "/assets/error_placeholder.png", [str(e)]
    finally:
        if tmp_path and os.path.exists(tmp_path):
            try:
                os.remove(tmp_path)
            except Exception:
                pass
