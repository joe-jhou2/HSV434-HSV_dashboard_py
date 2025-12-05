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

                # if (n_cell_types <= 1) {{
                #     width_factor = 2.5
                # }} else {{
                #     width_factor = 1.5
                # }}

                # if (n_genes <= 1) {{
                #     height_factor = 2.5
                # }} else {{
                #     height_factor = 1.5
                # }}

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


# utils/run_r_violin_plot.py

# import base64
# import io
# import time
# import warnings
# import matplotlib
# # Set backend to Agg for non-interactive server-side rendering
# matplotlib.use('Agg') 
# import matplotlib.pyplot as plt
# import seaborn as sns
# import pandas as pd
# import numpy as np

# def generate_violin_plot_from_df(df, color_file, selected_features, show_median_line=True):
#     """
#     V3: Optimized Python Violin Plot.
#     - Speed: Tuned KDE parameters (gridsize, cut) to fix 11s lag.
#     - Viz: Implements R-style 'facet_grid' with labels on Top and Left.
#     """
#     start_time = time.time()
#     print("--- Python Violin Plot (V3): Script execution started ---")
    
#     # Supress warnings
#     warnings.filterwarnings("ignore")
    
#     try:
#         # 1. PREPARE METADATA & COLORS
#         cell_type_colors = {k.replace(" ", "\n"): v for k, v in color_file.items()}
        
#         # 2. FILTER GENES
#         metadata_cols = {"Barcode", "UMAP_1", "UMAP_2", "CellType_Level3", "Subject", "Status"}
#         available_genes = list(set(df.columns) - metadata_cols)
#         target_genes = [g for g in selected_features if g in available_genes]

#         if not target_genes:
#             # Empty plot handler
#             fig, ax = plt.subplots(figsize=(6, 4))
#             ax.text(0.5, 0.5, "No data for filter selection", ha='center', va='center')
#             ax.axis('off')
#             buf = io.BytesIO()
#             plt.savefig(buf, format='png', dpi=150, bbox_inches='tight')
#             plt.close(fig)
#             buf.seek(0)
#             return f"data:image/png;base64,{base64.b64encode(buf.read()).decode()}", ["No data remaining."]

#         # 3. RESHAPE DATA
#         subset_cols = ["CellType_Level3", "Status"] + target_genes
#         plot_df = df[subset_cols].copy()
#         plot_df["CellType_Level3"] = plot_df["CellType_Level3"].astype(str).str.replace(" ", "\n")
        
#         melted_df = plot_df.melt(
#             id_vars=["CellType_Level3", "Status"], 
#             value_vars=target_genes,
#             var_name="Gene", 
#             value_name="Exp"
#         )
#         melted_df["Exp"] = pd.to_numeric(melted_df["Exp"], errors='coerce')

#         # Enforce Ordering
#         status_order = ["Prior", "Lesion", "Post"]
#         melted_df["Status"] = pd.Categorical(melted_df["Status"], categories=status_order, ordered=True)
#         melted_df["Gene"] = pd.Categorical(melted_df["Gene"], categories=target_genes, ordered=True)

#         # 4. DOWNSAMPLE (Keep this if dataset is massive, otherwise 10k is fine)
#         if len(melted_df) > 10000:
#             melted_df = melted_df.sample(n=10000, random_state=42)

#         # 5. GENERATE PLOT
#         n_genes = len(target_genes)
#         n_cell_types = melted_df["CellType_Level3"].nunique()
        
#         # Sizing
#         width_factor = 2.5 if n_cell_types <= 1 else 1.3
#         height_factor = 2.5 if n_genes <= 1 else 1.2
        
#         sns.set_theme(style="white", rc={"axes.facecolor": "white", "grid.color": "#f0f0f0"})
        
#         # --- KEY VIZ CHANGE: margin_titles=True ---
#         g = sns.FacetGrid(
#             melted_df, 
#             row="Gene", 
#             col="CellType_Level3", 
#             hue="CellType_Level3", 
#             palette=cell_type_colors,
#             sharex=True, 
#             sharey="row", 
#             height=height_factor, 
#             aspect=width_factor/height_factor,
#             margin_titles=True, # This puts labels on the edges, not inside every plot
#             despine=True
#         )

#         # --- KEY SPEED CHANGE: gridsize=30, cut=0 ---
#         # gridsize=30 (default 100) speeds up calculation 3x
#         # cut=0 prevents calculating tails beyond min/max data
#         g.map_dataframe(
#             sns.violinplot, 
#             x="Status", 
#             y="Exp", 
#             inner="quartile",
#             density_norm="width", 
#             linewidth=0.5,
#             alpha=0.8,
#             saturation=1,
#             gridsize=30, 
#             cut=0        
#         )

#         # Optional: Add Median Line (Fast version)
#         if show_median_line:
#             def fast_median_line(data, **kws):
#                 medians = data.groupby("Status", observed=False)["Exp"].median()
#                 plt.plot(medians.index, medians.values, color='black', linewidth=1, marker='.', markersize=3)
#             g.map_dataframe(fast_median_line)

#         # 6. FORMATTING (R-Style Labels)
        
#         # A. Clean Axis Labels
#         g.set_axis_labels("", "") 
        
#         # B. Handle Column Titles (Top)
#         g.set_titles(col_template="{col_name}", row_template="") # Clear row template first
        
#         # C. Handle Row Titles (Left Side Manually)
#         # Seaborn puts margin titles on the Right. We want Left.
#         # We iterate over the first column of axes to add the Gene label on the left.
#         genes_in_order = melted_df["Gene"].unique() # Should match row order
        
#         for ax, gene_name in zip(g.axes[:,0], genes_in_order):
#             # Add label to the left of the Y-axis
#             ax.annotate(gene_name, xy=(-0.3, 0.5), xycoords='axes fraction',
#                         ha='right', va='center', fontsize=11, fontweight='bold', rotation=0)
            
#             # Ensure Y labels (numbers) are visible only on left column
#             ax.tick_params(axis='y', labelleft=True)

#         # D. Ticks
#         for ax in g.axes.flat:
#             ax.tick_params(axis='x', rotation=45, labelsize=9)
#             ax.tick_params(axis='y', labelsize=8)
        
#         # Global Y Label
#         g.figure.text(0.02, 0.5, 'Log Normalized Expression', va='center', rotation='vertical', fontsize=12)
        
#         # Adjust layout to make room for the left labels we added
#         plt.subplots_adjust(left=0.2, bottom=0.15, right=0.98, top=0.92)

#         # 7. SAVE
#         buf = io.BytesIO()
#         g.savefig(buf, format='png', dpi=120, bbox_inches='tight') # dpi=120 is faster than 150
#         plt.close(g.fig)
#         buf.seek(0)
        
#         encoded_image = base64.b64encode(buf.read()).decode()
        
#         print(f"Violin Plot (V3) Computation time:  {round(time.time() - start_time, 3)}  sec")
#         return f"data:image/png;base64,{encoded_image}", []

#     except Exception as e:
#         print(f"--- ERROR generating Violin Plot (Python) ---\n{e}")
#         return "/assets/error_placeholder.png", [str(e)]