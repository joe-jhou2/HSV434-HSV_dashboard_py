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
                    cat("⚠️ selectedCellType empty — defaulting to 'All'\n", file = stderr())
                    selectedCellType <- "All"
                }}
                if (length(selectedFeature) == 0) {{
                    stop("❌ No selected features provided from Python.")
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

                    # if (n_cell_types <= 1) {{
                    #     height_factor = 2.5
                    # }} else {{
                    #     height_factor = 0.3
                    # }}

                    # if (n_genes <= 1) {{
                    #     width_factor = 1
                    # }} else {{
                    #     width_factor = 0.65
                    # }}

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

# utils/run_r_dot_plot.py

# import base64
# import io
# import time
# import matplotlib
# # Set backend to Agg for non-interactive server-side rendering
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# import matplotlib.colors as mcolors
# import seaborn as sns
# import pandas as pd
# import numpy as np

# def generate_dot_plot_from_df(data_pert, data_gex, color_file, selected_features, selected_celltypes):
#     """
#     Native Python implementation of the Dot Plot.
#     - Replaces R computation with Pandas (Vectorized).
#     - Replaces ggplot2 with Seaborn/Matplotlib.
#     - Fixes the 'object not found' crash by calculating dimensions locally.
#     """
#     start_time = time.time()
#     print("--- Python Dot Plot: Script execution started ---")

#     try:
#         # 1. VALIDATION & SETUP
#         if not selected_features:
#             raise ValueError("No selected features provided.")
        
#         # Handle 'All' selection for cell types
#         # (If selected_celltypes is empty or contains "All", we use all available in data)
#         available_celltypes = data_pert['CellType_Level3'].unique()
#         if not selected_celltypes or "All" in selected_celltypes:
#              target_celltypes = available_celltypes
#         else:
#              target_celltypes = [ct for ct in selected_celltypes if ct in available_celltypes]

#         if len(target_celltypes) == 0:
#             raise ValueError("No matching cell types found in data.")

#         # 2. DATA WRANGLING (Pandas > R Dplyr)
#         # ---------------------------------------------------------
        
#         # A. Process Percent Data (Size of Dot)
#         # -------------------------------------
#         # Filter cols: Subject, CellType, Status + Genes
#         pert_cols = ["Subject", "CellType_Level3", "Status"] + [f for f in selected_features if f in data_pert.columns]
#         df_pert = data_pert[pert_cols].copy()
#         df_pert = df_pert.fillna(0) # Match R: GenePert[is.na(GenePert)] = 0
        
#         # Melt to Long Format
#         pert_long = df_pert.melt(
#             id_vars=["Subject", "CellType_Level3", "Status"],
#             var_name="Gene", value_name="percent"
#         )
        
#         # Filter by CellType
#         pert_long = pert_long[pert_long["CellType_Level3"].isin(target_celltypes)]
        
#         # Aggregate: Mean Percent per (CellType, Status, Gene)
#         agg_pert = pert_long.groupby(["CellType_Level3", "Status", "Gene"], observed=False)["percent"].mean().reset_index()
#         agg_pert.rename(columns={"percent": "Avg_pert"}, inplace=True)

#         # B. Process Expression Data (Color of Dot)
#         # -----------------------------------------
#         # Filter cols
#         gex_cols = ["Subject", "CellType_Level3", "Status"] + [f for f in selected_features if f in data_gex.columns]
#         df_gex = data_gex[gex_cols].copy()
        
#         # R Logic: Group by (Subject_Status, CellType) -> Mean -> Then separate Subject_Status
#         # In Python, we can just GroupBy [Subject, Status, CellType] directly
#         # Calculate mean expression per Subject first (Sample Level)
#         subject_means = df_gex.groupby(["Subject", "Status", "CellType_Level3"], observed=False).mean().reset_index()
        
#         # Melt Subject Means to Long
#         gex_long = subject_means.melt(
#             id_vars=["Subject", "Status", "CellType_Level3"],
#             var_name="Gene", value_name="expression"
#         )
        
#         # Filter by CellType
#         gex_long = gex_long[gex_long["CellType_Level3"].isin(target_celltypes)]
        
#         # Aggregate: Mean Expression per (CellType, Status, Gene) across subjects
#         agg_gex = gex_long.groupby(["CellType_Level3", "Status", "Gene"], observed=False)["expression"].mean().reset_index()
#         agg_gex.rename(columns={"expression": "Avg_exp"}, inplace=True)

#         # 3. MERGE & ORDERING
#         # ---------------------------------------------------------
#         plot_data = pd.merge(agg_pert, agg_gex, on=["CellType_Level3", "Status", "Gene"], how="inner")

#         # Enforce Orders
#         # Status
#         status_order = ["Prior", "Lesion", "Post"]
#         plot_data["Status"] = pd.Categorical(plot_data["Status"], categories=status_order, ordered=True)
        
#         # Genes (Match user selection order)
#         plot_data["Gene"] = pd.Categorical(plot_data["Gene"], categories=selected_features, ordered=True)
        
#         # Cell Types (Match color file keys order)
#         # Clean keys (spaces -> newlines)
#         cleaned_color_keys = [k.replace(" ", "\n") for k in color_file.keys()]
#         plot_data["CellType_Level3"] = plot_data["CellType_Level3"].astype(str).str.replace(" ", "\n")
        
#         # Filter ordered keys to only those present in data to avoid empty rows
#         present_types = plot_data["CellType_Level3"].unique()
#         final_type_order = [ct for ct in cleaned_color_keys if ct in present_types]
        
#         # If any types are missing from color keys, append them
#         remaining = [ct for ct in present_types if ct not in final_type_order]
#         final_type_order.extend(remaining)
        
#         plot_data["CellType_Level3"] = pd.Categorical(
#             plot_data["CellType_Level3"], 
#             categories=final_type_order, # Reverse to plot top-down if using Y-axis
#             ordered=True
#         )

#         if plot_data.empty:
#             raise ValueError("No data available for the selected features and cell types.")

#         # 4. PLOTTING
#         # ---------------------------------------------------------
#         n_genes = plot_data["Gene"].nunique()
#         n_cell_types = plot_data["CellType_Level3"].nunique()

#         # Dynamic Sizing (Matching R Logic roughly)
#         # R used: width = n_genes * 0.65, height = n_cell_types * 0.3
#         # Matplotlib needs slightly larger base inches
#         fig_width = max(n_genes * 0.8 + 2, 4) 
#         fig_height = max(n_cell_types * 0.4 + 1, 3)

#         sns.set_theme(style="white", rc={"axes.facecolor": "white", "grid.color": "#ededed"})
        
#         # Define Custom Colormap (Gray -> Red)
#         cmap = mcolors.LinearSegmentedColormap.from_list("gray_red", ["#E0E0E0", "red"])

#         # Create Grid: Columns = Genes
#         g = sns.FacetGrid(
#             plot_data, 
#             col="Gene", 
#             sharex=True, 
#             sharey=True,
#             height=fig_height,
#             aspect=(fig_width / n_genes) / fig_height if n_genes > 0 else 1
#         )

#         # Draw Scatter Plots
#         # size_norm sets the scaling for the bubble sizes (0 to 100 percent)
#         g.map_dataframe(
#             sns.scatterplot,
#             x="Status",
#             y="CellType_Level3",
#             size="Avg_pert",
#             hue="Avg_exp",
#             sizes=(10, 200), # Min/Max pixel size of dots
#             size_norm=(0, 100),
#             palette=cmap,
#             edgecolor="black",
#             linewidth=0.5,
#             alpha=1
#         )

#         # 5. FORMATTING & LEGEND
#         # ---------------------------------------------------------
#         g.set_axis_labels("", "") # Hide X/Y labels
#         g.set_titles(col_template="{col_name}") # Clean Gene Titles
        
#         # Ticks formatting
#         for ax in g.axes.flat:
#             ax.tick_params(axis='x', rotation=90, labelsize=9)
#             ax.tick_params(axis='y', labelsize=9)
#             # Add light grid for readability
#             ax.grid(True, axis='y', linestyle='--', alpha=0.5)

#         # Adjust Layout
#         plt.subplots_adjust(top=0.85, bottom=0.2, left=0.15, right=0.95)

#         # Manually Create Legends (Seaborn FacetGrid legend is often tricky)
#         # We assume the last axis has the plot properties
#         # Remove default legend
#         g.fig.legends = []
        
#         # 1. Color Bar (Avg Expression)
#         cbar_ax = g.fig.add_axes([0.96, 0.3, 0.015, 0.4]) # [left, bottom, width, height]
#         norm = plt.Normalize(plot_data['Avg_exp'].min(), plot_data['Avg_exp'].max())
#         cb = plt.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cbar_ax)
#         cb.set_label("Avg Expression", size=8)
#         cb.ax.tick_params(labelsize=8)

#         # 2. Size Legend (Percent Expressed) - approximate manually or use a dummy plot
#         # For simplicity in FacetGrid, we let the user hover or interpret relative size
#         # Or we can add a simple text note
#         g.fig.text(0.96, 0.25, "Size: % Exp", ha="center", fontsize=8)

#         # 6. SAVE
#         buf = io.BytesIO()
#         g.savefig(buf, format='png', dpi=150, bbox_inches='tight')
#         plt.close(g.fig)
#         buf.seek(0)
#         encoded_image = base64.b64encode(buf.read()).decode()

#         duration = round(time.time() - start_time, 3)
#         print(f"Dot Plot Computation time:  {duration}  sec")
        
#         return f"data:image/png;base64,{encoded_image}", []

#     except Exception as e:
#         print(f"--- ERROR generating Dot Plot (Python) ---\n{e}")
#         return "/assets/error_placeholder.png", [str(e)]