# pages/tabs/scrnaseq_gene_tab.py

import os
from dash import html, dcc, Input, Output, State
import dash
from dotenv import load_dotenv
import dash_bootstrap_components as dbc
from utils.helper import get_dataset_options
from utils.gene_utils import (set_refresh_flag, get_refresh_flag,
                             get_available_gene_universe,
                             get_gene_list,
                             check_genes_availability,
                             build_ordered_gene_list,
                             run_precompute_r_async,
                             log_progress)
from utils.run_gex_data_loader import load_filtered_gex_data
from utils.run_pert_data_loader import load_filtered_pert_data
from utils.run_r_gene_umap import generate_feature_umap_from_df
from utils.run_r_gene_heatmap import generate_heatmap_from_df
from utils.run_r_gene_violin_plot import generate_violin_plot_from_df
from utils.run_r_gene_dot_plot import generate_dot_plot_from_df

# Load environment variables from .env file
load_dotenv()

# --- Layout Definition ---
gene_tab_layout = html.Div([
    html.H4("Gene Discovery Controls"),

    dbc.Row([
        dbc.Col([
            html.Label("Select Clusters"),
            dcc.Dropdown(id="cluster-dropdown-gene-tab", placeholder="Default: All", multi=True)
        ], width=4),
        dbc.Col([
            html.Label("Select Subjects"),
            dcc.Dropdown(id="subject-dropdown-gene-tab", placeholder="Default: All", multi=True)
        ], width=4),
        dbc.Col([
            html.Label("Select Genes (or type)"),
            dcc.Dropdown(id="gene-dropdown-gene-tab", multi=True, placeholder="Select from list ..."),
            dcc.Input(id="gene-input-gene-tab", type="text", placeholder="Or type gene names separated by commas", style={"marginTop": "10px"})
        ], width=4)
    ]),

    dbc.Button("Run Analysis", id="run-gene-tab-btn", color="primary", className="my-3"),
    html.Hr(),

    # --- Add Sub-Tabs ---
    dcc.Tabs(id="gene-sub-tabs", value="feature-umap-sub-tab", children=[

        # 1. Feature UMAP Sub-Tab
        dcc.Tab(label="Feature UMAP", value="feature-umap-sub-tab", className="feature-umap-sub-tab", children=[
            html.Div(className="py-4", children=[
                dbc.Alert(id="feature-umap-notify", color="warning", is_open=False, duration=4000),
                dcc.Loading(type="circle", children=[
                    dbc.Row([
                        dbc.Col(html.Img(id="feature-umap-plot-img", style={'width': '100%'}))
                    ])
                ])
            ])
        ]),

        # 2. Heatmap Sub-Tab
        dcc.Tab(label="Heatmap", value="heatmap-sub-tab", className="heatmap-sub-tab", children=[
            html.Div(className="py-4", children=[
                dcc.Loading(type="circle", children=[
                    dbc.Row([
                        dbc.Col(html.Img(id="heatmap-plot-img", style={'width': '100%'}))
                    ])
                ])
            ])
        ]),

        # 3. Violin Plots Sub-Tab
        dcc.Tab(label="Violin Plots", value="vln-sub-tab", className="vln-sub-tab", children=[
             html.Div(className="py-4", children=[
                 dcc.Loading(type="circle", children=[
                    dbc.Row([
                        dbc.Col(html.Img(id="violin-plots-img", style={'width': '100%'}))
                    ])
                 ])
             ])
        ]),

        # 4. Dot Plots Sub-Tab
        dcc.Tab(label="Dot Plots", value="dot-sub-tab", className="dot-sub-tab", children=[
             html.Div(className="py-4", children=[
                 dcc.Loading(type="circle", children=[
                    dbc.Row([
                        dbc.Col(html.Img(id="dot-plots-img", style={'width': '100%'}))
                    ])
                 ])
             ])
        ]),
    ]),

    dcc.Interval(
    id="gene-refresh-interval",
    interval=10 * 1000,  # every 10 seconds
    n_intervals=0
)
])

# --- Callback Registration ---
def register_callbacks(app):
    # Callback 1: Update all dropdown options when the main dataset selection changes
    @app.callback(
        Output("cluster-dropdown-gene-tab", "options"),
        Output("subject-dropdown-gene-tab", "options"),
        Output("gene-dropdown-gene-tab", "options"),
        Output("cluster-dropdown-gene-tab", "value"),
        Output("subject-dropdown-gene-tab", "value"),
        Output("gene-dropdown-gene-tab", "value"),
        Input("dataset_option", "value")
    )
    def update_gene_tab_dropdowns(dataset_prefix):
        if not dataset_prefix:
            return [], [], [], [], [], []
        
        options = get_dataset_options(dataset_prefix, bucket_name=None, force_s3=True)
        gene_list = get_gene_list(dataset_prefix, bucket_name=None, force_s3=True)
        
        cluster_options = ['All'] + options["clusters"]
        subject_options = ['All'] + options["subjects"]
        
        # Set some default genes for the user
        default_genes = ['CD4', 'CD8A', 'CD14', 'FCGR3A', 'IFNG', 'PRF1','GZMA', 'GZMB']
        # Ensure default genes actually exist in the list
        valid_default_genes = [g for g in default_genes if g in gene_list]

        return cluster_options, subject_options, gene_list, ['All'], ['All'], valid_default_genes

    # Callback 2: Generate plots when the "Run Analysis" button is clicked
    @app.callback(
        Output("feature-umap-plot-img", "src", allow_duplicate=True),
        Output("heatmap-plot-img", "src", allow_duplicate=True), 
        Output("violin-plots-img", "src", allow_duplicate=True),
        Output("dot-plots-img", "src", allow_duplicate=True),
        Output("feature-umap-notify", "children", allow_duplicate=True),  # Output for notification
        Output("feature-umap-notify", "is_open", allow_duplicate=True),   # Output to show/hide notification
        Input("run-gene-tab-btn", "n_clicks"),
        State("dataset_option", "value"),
        State("cluster-dropdown-gene-tab", "value"),
        State("subject-dropdown-gene-tab", "value"),
        State("gene-dropdown-gene-tab", "value"),
        State("gene-input-gene-tab", "value"),
        prevent_initial_call=True
    )
    def update_gene_tab_plots(n_clicks, dataset_prefix, selected_clusters, selected_subjects, selected_genes, typed_genes):
        if not dataset_prefix:
            return "", "", "", "", "Please select a dataset and genes to visualize.", True

        # --- Combine dropdown and typed input ---
        requested_genes = build_ordered_gene_list(selected_genes, typed_genes)

        if not requested_genes:
            return "", "", "", "", "Please select or type genes to visualize.", True

        # --- Validate against master available gene list ---
        valid_gene_universe = set(get_available_gene_universe(dataset_prefix, bucket_name=None, force_s3=True))
        if not valid_gene_universe:
            return "", "", "", "", (
                f"No available gene list found for {dataset_prefix}. "
                "Please generate it first using export_available_genes.R."
            ), True

        invalid_genes = [g for g in requested_genes if g not in valid_gene_universe]
        valid_genes = [g for g in requested_genes if g in valid_gene_universe]

        if invalid_genes:
            msg_invalid = f"Unrecognized genes ignored: {', '.join(invalid_genes)}"
            return "", "", "", "", f"{msg_invalid}\n. Please check input.", True
        else:
            msg_invalid = ""

        if not valid_genes:
            return "", "", "", "", f"{msg_invalid}\n No valid genes to visualize. Please check input.", True

        # Continue downstream only with valid genes
        requested_genes = valid_genes
        
        # Normalize selections
        clusters_to_filter = [] if not selected_clusters or "All" in selected_clusters else selected_clusters
        subjects_to_filter = [] if not selected_subjects or "All" in selected_subjects else selected_subjects

        # --- Step 1: Check local JSON for available genes ---
        genes_available, missing_genes = check_genes_availability(dataset_prefix, requested_genes, bucket_name=None, force_s3=True)
        notification_msg = ""
        notification_open = False
        
        # --- Step 2: Handle missing genes ---
        if missing_genes:
            log_progress(f"Genes missing: {missing_genes}")
            set_refresh_flag(dataset_prefix, "running")
            run_precompute_r_async(dataset_prefix, missing_genes)

            if not genes_available:
                log_progress("No local genes available; skipping immediate plot generation.")
    
                # Show a simple placeholder
                placeholder_path = "/assets/images/HSV.png"
                placeholder_src = f"/{placeholder_path}" if os.path.exists(placeholder_path) else ""

                notification_msg = (
                    f"Retrieving {len(missing_genes)} missing genes "
                    f"({', '.join(missing_genes)}) from DataLake... "
                    "Plots will refresh automatically when ready."
                )
                notification_open = True

                # Start async background pull
                set_refresh_flag(dataset_prefix, "running")
                run_precompute_r_async(dataset_prefix, missing_genes)

                log_progress("Background retrieval started; UI stays interactive.")
                return placeholder_src, placeholder_src, placeholder_src, placeholder_src, notification_msg, notification_open

            # Partial case (some missing, some available)
            notification_msg = (
                    f"Fetching {len(missing_genes)} missing gene(s): {', '.join(missing_genes)}. "
                    "Plots will auto-refresh once retrieval completes."
                )
            notification_open = True

        # --- Step 3: Proceed to plotting with available genes ---
        umap_src = ""
        heatmap_src = ""
        violin_src = ""
        dot_src = ""
        
        if genes_available:
            try:
                # Query:
                data_gex, color_map = load_filtered_gex_data(
                    dataset_prefix,
                    genes=genes_available,
                    clusters=clusters_to_filter,
                    subjects=subjects_to_filter,
                    bucket_name=None,
                    force_s3=True
                )

                data_pert, color_map = load_filtered_pert_data(
                    dataset_prefix,
                    genes=genes_available,
                    clusters=clusters_to_filter,
                    subjects=subjects_to_filter,
                    bucket_name=None,
                    force_s3=True
                )

                umap_src, _ = generate_feature_umap_from_df(data_gex, genes_available)
                heatmap_src, _ = generate_heatmap_from_df(data_gex, color_map, genes_available)
                violin_src, _ = generate_violin_plot_from_df(data_gex, color_map, genes_available)
                dot_src, _ = generate_dot_plot_from_df(data_pert, data_gex, color_map, genes_available, clusters_to_filter)
                
            except Exception as e:
                log_progress(f"Error generating plots: {e}")
                umap_src = "/assets/images/HSV.png"
                heatmap_src = "/assets/images/HSV.png"
                violin_src = "/assets/images/HSV.png"
                dot_src = "/assets/images/HSV.png"
        else:
            # Only placeholder if no data locally
            umap_src = "/assets/images/HSV.png"
            heatmap_src = "/assets/images/HSV.png"
            violin_src = "/assets/images/HSV.png"
            dot_src = "/assets/images/HSV.png"

        return umap_src, heatmap_src, violin_src, dot_src, notification_msg, notification_open

    # # Callback 3: Periodic refresh to check for completed background jobs
    @app.callback(
        Output("feature-umap-plot-img", "src", allow_duplicate=True),
        Output("heatmap-plot-img", "src", allow_duplicate=True),
        Output("violin-plots-img", "src", allow_duplicate=True),
        Output("dot-plots-img", "src", allow_duplicate=True),
        Output("feature-umap-notify", "children", allow_duplicate=True),
        Output("feature-umap-notify", "is_open", allow_duplicate=True),
        Input("gene-refresh-interval", "n_intervals"),
        State("dataset_option", "value"),
        State("cluster-dropdown-gene-tab", "value"),
        State("subject-dropdown-gene-tab", "value"),
        State("gene-dropdown-gene-tab", "value"),
        State("gene-input-gene-tab", "value"),
        prevent_initial_call=True
    )

    def auto_refresh_gene_data(n, dataset_prefix, selected_clusters, selected_subjects,
                           selected_genes, typed_genes):

        if not dataset_prefix:
            raise dash.exceptions.PreventUpdate

        # Only refresh when R script completed successfully
        if get_refresh_flag(dataset_prefix) != "ready":
            raise dash.exceptions.PreventUpdate

        log_progress(f"Auto-refresh triggered for {dataset_prefix}")

        # rebuild gene list fresh, not from stale callback state
        requested_genes = build_ordered_gene_list(selected_genes, typed_genes)

        # re-check availability
        genes_available, missing = check_genes_availability(dataset_prefix, requested_genes, bucket_name=None, force_s3=True)

        # now ALL missing genes should exist → safest to use requested_genes directly
        final_genes = genes_available or requested_genes

        clusters_to_filter = [] if not selected_clusters or "All" in selected_clusters else selected_clusters
        subjects_to_filter = [] if not selected_subjects or "All" in selected_subjects else selected_subjects

        try:
            data_gex, color_map = load_filtered_gex_data(
                dataset_prefix,
                genes=final_genes,
                clusters=clusters_to_filter,
                subjects=subjects_to_filter,
                bucket_name=None,
                force_s3=True
            )

            data_pert, color_map = load_filtered_pert_data(
                dataset_prefix,
                genes=final_genes,
                clusters=clusters_to_filter,
                subjects=subjects_to_filter,
                bucket_name=None,
                force_s3=True
            )

            umap_src, _ = generate_feature_umap_from_df(data_gex, final_genes)
            heatmap_src, _ = generate_heatmap_from_df(data_gex, color_map, final_genes)
            violin_src, _ = generate_violin_plot_from_df(data_gex, color_map, final_genes)
            dot_src, _ = generate_dot_plot_from_df(data_pert, data_gex, color_map, final_genes, clusters_to_filter)

            msg = f"New genes added for {dataset_prefix}. Plots updated."

        except Exception as e:
            log_progress(f"Auto-refresh error: {e}")
            return "", "", "", "", f"Auto-refresh failed: {e}", True

        # After successful refresh, reset flag
        set_refresh_flag(dataset_prefix, "idle")

        return umap_src, heatmap_src, violin_src, dot_src, msg, True
