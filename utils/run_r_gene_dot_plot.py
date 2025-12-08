# utils/run_r_dot_plot.py

import base64
import os
import tempfile
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from utils.db_connection import dict_to_r_vector

def generate_dot_plot_from_df(data_pert, data_gex, color_file, selected_features, selected_celltypes):
    """
    Generate a dot plot using the R computeDotPlot() function from a pandas DataFrame.
    """
    tmp_path = ""
    try:
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            tmp_path = tmp.name
        tmp_path_r = tmp_path.replace("\\", "/")

        r_color_string = dict_to_r_vector(color_file)

        # Pass the dataframe to R
        with localconverter(ro.default_converter + pandas2ri.converter):
            ro.globalenv["data_pert_py"] = data_pert
            ro.globalenv["data_gex_py"] = data_gex
            # ro.globalenv["color_file"] = color_file
            ro.globalenv['selected_features_r'] = ro.StrVector(selected_features)
            ro.globalenv['selected_celltypes_r'] = ro.StrVector(selected_celltypes)
        
        # color_file_r = color_file.replace("\\", "/")

        r_code = f"""
            suppressPackageStartupMessages({{
                library(ggplot2)
                library(dplyr)
                library(reshape2)
                library(jsonlite)
                library(tidyr)
                library(ggh4x)
            }})

            # Suppressing warnings speeds up execution by preventing console I/O lag
            suppressWarnings(suppressMessages({{
            
                # Start timing
                start_time <- Sys.time()
                cat("--- R Dot Plot: Script execution started ---\\n", file=stderr())

                CellType_color <- {r_color_string}
                if (is.list(CellType_color)) {{
                    CellType_color <- unlist(CellType_color)
                }}
    
                # --- Assign Python variables in R ---
                data_pert = data_pert_py
                data_gex = data_gex_py
                selectedFeature = selected_features_r
                selectedCellType = selected_celltypes_r

                if (length(selectedCellType) == 0) {{
                    cat("selectedCellType empty — defaulting to 'All'\n", file = stderr())
                    selectedCellType <- "All"
                }}
                if (length(selectedFeature) == 0) {{
                    stop("No selected features provided from Python.")
                }}

                suppressWarnings({{
                    # Calculate Cytokine+ Pert within CellType and Status
                    GenePert = data_pert
                    GenePert[is.na(GenePert)] = 0
                    
                    # Extracte Cytokine/Gene expression data
                    Pert_Gene = GenePert %>%
                        dplyr::select(Subject, CellType_Level3, Status, all_of(selectedFeature)) %>% 
                        group_by(CellType_Level3, Status) %>%
                        dplyr::filter(if (any(selectedCellType == "All")) TRUE
                                        else CellType_Level3 %in% selectedCellType) %>%
                        droplevels() %>%
                        pivot_longer(cols = all_of(selectedFeature), names_to = "Gene", values_to = "percent") %>%
                        dplyr::mutate(Gene = factor(Gene, levels = selectedFeature))
                    
                    # Calculate Gene Intensity
                    Exp_Gene = data_gex %>%
                        dplyr::select(Subject, CellType_Level3, Status, all_of(selectedFeature)) %>%
                        mutate(orig.ident = paste(Subject, Status, sep = "_")) %>%
                        group_by(orig.ident, CellType_Level3) %>%
                        summarise(across(-c(Subject, Status), function(x) mean(x, na.rm = TRUE))) %>%
                        separate(orig.ident, into = c("Subject", "Status"), sep = "_") %>%
                        dplyr::filter(if (any(selectedCellType == "All")) TRUE
                                        else CellType_Level3 %in% selectedCellType) %>%
                        droplevels() %>% 
                        pivot_longer(cols = all_of(selectedFeature), names_to = "Gene", values_to = "expression") %>%
                        dplyr::mutate(Gene = factor(Gene, levels = selectedFeature))
                        
                    # Prepare data for plotting
                    data.plot1 = Pert_Gene %>% 
                        group_by(CellType_Level3, Status, Gene) %>%
                        summarise(Avg_pert = mean(percent), .groups = 'drop') %>%
                        dplyr::mutate(idx = paste(CellType_Level3, Status, Gene, sep = "_"))

                    data.plot2 = Exp_Gene %>%
                        group_by(CellType_Level3, Status, Gene) %>%
                        summarise(Avg_exp = mean(expression), .groups = 'drop') %>%
                        dplyr::mutate(idx = paste(CellType_Level3, Status, Gene, sep = "_"))

                    plot.data = data.plot1 %>% left_join(data.plot2[,c("idx", "Avg_exp")], by = "idx") %>%
                        dplyr::mutate(Gene = factor(Gene, levels = selectedFeature),
                                    Status = factor(Status, levels = c("Prior", "Lesion", "Post")),
                                    CellType_Level3 = factor(CellType_Level3, levels = names(CellType_color)))
                    
                    # plotting
                    if (nrow(plot.data) == 0) {{
                        stop("No data available for the selected features and cell types.")
                    }} else {{
                        g = ggplot(plot.data, aes(x = Status, y = CellType_Level3)) + 
                            geom_point(aes(fill = Avg_exp, size = Avg_pert), colour = "black", pch = 21, stroke = 0.1) +
                            scale_fill_gradient2(mid = "gray", high = "red") +
                            facet_grid2(CellType_Level3 ~ Gene, scales = "free_y") +
                            theme_minimal(base_size = 8) + 
                            theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.2),
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(),
                                    panel.background = element_blank(),
                                    strip.text = element_text(size = 8),
                                    strip.text.x = element_text(size = 8, colour = "black", angle = 0, margin = margin(t = 0, r = 0, b = 2, l = 0)),
                                    strip.text.y = element_text(size = 0, colour = "black", angle = 0, margin = margin(t = 0, r = 0, b = 0, l = 0)),
                                    axis.ticks.length = unit(0, "pt"),
                                    axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5, colour = "black", margin = margin(t = 2, r = 0, b = 0, l = 0)),
                                    axis.text.y = element_text(size = 8, colour = "black", margin = margin(t = 0, r = 2, b = 0, l = 0)),
                                    axis.title.x = element_blank(),
                                    axis.title.y = element_blank(),
                                    plot.margin = margin(0,0,0,0),
                                    legend.box = "vertical",
                                    legend.position = "bottom",
                                    legend.box.just = "center",
                                    legend.margin = margin(0, 0, 0, 0),
                                    legend.box.margin = margin(0, 0, 0, 0),
                                    legend.title = element_text(size = 8),
                                    legend.text = element_text(size = 8),
                                    legend.key.height = unit(1.5, "lines"),
                                    legend.key.width = unit(1.5, "lines"),
                                    panel.spacing = unit(0.02, "lines"))  +
                            scale_size_continuous(range = c(0.1, 5), limits = c(0, 100)) + 
                            labs(size = "Percent Expressed", fill = "Average Expression") +
                            guides(
                                    size = guide_legend(
                                        title.position = "top",
                                        title.hjust = 0.5,
                                        nrow = 1, 
                                        override.aes = list(shape = 21, colour = "black", fill = "gray", stroke = 0.2),
                                        keyheight = unit(1, "lines"),
                                        keywidth  = unit(1.5, "lines")
                                    ),
                                    fill = guide_colorbar(
                                        title.position = "top",
                                        title.hjust = 0.5
                                    )
                                )
                    }}

                    # --- Save plot ---
                    n_genes = length(unique(exp_meta_df_reshape$Gene))
                    n_cell_types = length(unique(exp_meta_df_reshape$CellType_Level3))

                    base_w = 4
                    base_h = 3
                    
                    # Add incremental width/height
                    calc_w = base_w + (n_genes * 0.7)
                    calc_h = base_h + (n_cell_types * 0.3)

                    ggsave("{tmp_path_r}", plot = g, device = "png", width = calc_w, height = calc_h, dpi = 300, limitsize = FALSE)

                    # End timing
                    end_time <- Sys.time()
                    cat("DotPlot Computation time: ", round(end_time - save_start, 3), " sec\n")
                }})
            }}))
        """
        ro.r(r_code)

        # Encode and return
        with open(tmp_path, "rb") as image_file:
            encoded_image = base64.b64encode(image_file.read()).decode()
        return f"data:image/png;base64,{encoded_image}", []

    except Exception as e:
        print(f"--- ERROR generating Dot Plot from df ---\\n{e}")
        return "/assets/error_placeholder.png", [str(e)]
    finally:
        if tmp_path and os.path.exists(tmp_path):
            try:
                os.remove(tmp_path)
            except Exception:
                pass