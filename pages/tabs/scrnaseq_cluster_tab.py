# pages/tabs/scrnaseq_cluster_tab.py

from dash import html, dcc, Input, Output, State
from dotenv import load_dotenv
import dash_bootstrap_components as dbc
from utils.helper import get_dataset_options
from utils.run_r_cluster_umap import generate_umap_plot
from utils.run_r_cluster_stat import generate_clusterStat_plots
from utils.run_r_cluster_perSubject import generate_PerSubject_StackBar_plots

# Load environment variables from .env file
load_dotenv()

cluster_tab_layout = html.Div([
    html.H4("Cluster Discovery Controls"),
    dbc.Row([
        dbc.Col([
            html.Label("Select Clusters (CellType_Level3)"),
            dcc.Dropdown(id="cluster-dropdown-cluster-tab", placeholder="Default: All", multi=True)
        ], width=4),
        dbc.Col([
            html.Label("Select Subjects"),
            dcc.Dropdown(id="subject-dropdown-cluster-tab", placeholder="Default: All", multi=True)
        ], width=4)
    ]),
    dbc.Button("Run Analysis", id="run-cluster-btn", color="primary", className="my-3"),
    html.Hr(),

    # --- Here is the new nested Tabs component ---
    dcc.Tabs(id="cluster-sub-tabs", value="umap-sub-tab", children=[
        
        # 1. UMAP Plots Sub-Tab
        dcc.Tab(label="UMAP Plots", value="umap-sub-tab", className="umap-sub-tab", children=[
            html.Div(className="py-4", children=[
                dcc.Loading(type="circle", children=[
                    dbc.Row([
                        dbc.Col(html.Img(id="cluster-umap-all-img", className="img-fluid"), width=6, className="mx-auto")
                    ]),
                    dbc.Row([
                        dbc.Col(html.Img(id="cluster-umap-prior-img", className="img-fluid"), width=4),
                        dbc.Col(html.Img(id="cluster-umap-lesion-img", className="img-fluid"), width=4),
                        dbc.Col(html.Img(id="cluster-umap-post-img", className="img-fluid"), width=4)
                    ], className="mt-3 g-3")
                ])
            ])
        ]),
        
        # 2. Cell Type Stats Sub-Tab
        dcc.Tab(label="Cell Type Stats", value="stats-sub-tab", className="stats-sub-tab", children=[
            html.Div(className="py-4", children=[
                dcc.Loading(type="circle", children=[
                    dbc.Row([
                        dbc.Col(html.Img(id="celltype-stats-barplot-img", style={'width': '100%'}), width=12),
                    ], className="g-3")
                ])
            ])
        ]),
        
        # 3. Cell Type by Subject Sub-Tab
        dcc.Tab(label="Cell Type by Subject", value="subject-sub-tab", className="subject-sub-tab", children=[
            html.Div(className="py-4", children=[
                 dcc.Loading(type="circle", children=[
                    dbc.Row([
                        dbc.Col(html.Img(id="celltype-subject-barplot-img", style={'width': '100%'}), width=12),
                    ])
                ])
            ])
        ]),
    ])
])

# --- Final Callback with State and Plotting Logic ---
def register_callbacks(app):
    # This first callback updates the dropdown options when the dataset changes
    @app.callback(
        Output("cluster-dropdown-cluster-tab", "options"),
        Output("subject-dropdown-cluster-tab", "options"),
        Output("cluster-dropdown-cluster-tab", "value"),
        Output("subject-dropdown-cluster-tab", "value"),
        Input("dataset_option", "value")
    )
    def update_dropdown_options(dataset_prefix):
        if not dataset_prefix:
            return [], [], [], []
        options = get_dataset_options(dataset_prefix, bucket_name=None, force_s3=True)
        
        # Prepend 'All' to the lists of options
        cluster_options = ['All'] + options["clusters"]
        subject_options = ['All'] + options["subjects"]

        return cluster_options, subject_options, ['All'], ['All']

    # This second callback generates the plots
    @app.callback(
        Output("cluster-umap-all-img", "src"),
        Output("cluster-umap-prior-img", "src"),
        Output("cluster-umap-lesion-img", "src"),
        Output("cluster-umap-post-img", "src"),

        Output("celltype-stats-barplot-img", "src"),

        Output("celltype-subject-barplot-img", "src"),
        
        Input("run-cluster-btn", "n_clicks"),
        State("dataset_option", "value"),
        State("cluster-dropdown-cluster-tab", "value"),
        State("subject-dropdown-cluster-tab", "value"),
        prevent_initial_call=True
    )
    def update_cluster_tab_plots(n_clicks, dataset_prefix, selected_clusters, selected_subjects):
        if not dataset_prefix:
            return [""] * 4
        
        clusters_to_filter = [] if not selected_clusters or 'All' in selected_clusters else selected_clusters
        subjects_to_filter = [] if not selected_subjects or 'All' in selected_subjects else selected_subjects
        
        # Pass the selected clusters and subjects to the plotting function
        src_all = generate_umap_plot(dataset_prefix, "All", "All Timepoints", clusters_to_filter, subjects_to_filter)
        src_prior = generate_umap_plot(dataset_prefix, "Prior", "Prior", clusters_to_filter, subjects_to_filter)
        src_lesion = generate_umap_plot(dataset_prefix, "Lesion", "Lesion", clusters_to_filter, subjects_to_filter)
        src_post = generate_umap_plot(dataset_prefix, "Post", "Post", clusters_to_filter, subjects_to_filter)
        
        stats_plot_src = generate_clusterStat_plots(dataset_prefix)

        subject_plot_src = generate_PerSubject_StackBar_plots(dataset_prefix, subjects=subjects_to_filter) # Pass the filtered list

        return src_all, src_prior, src_lesion, src_post, stats_plot_src, subject_plot_src