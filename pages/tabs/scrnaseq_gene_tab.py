# pages/tabs/scrnaseq_gene_tab.py

import os
import subprocess
import datetime
import time
from dash import html, dcc, Input, Output, State
import dash
import pandas as pd
import json
from io import BytesIO
import boto3
from urllib.parse import urlparse
from dotenv import load_dotenv
import dash_bootstrap_components as dbc
from utils.run_gex_data_loader import load_filtered_gex_data
from utils.run_pert_data_loader import load_filtered_pert_data
from utils.run_r_gene_umap import generate_feature_umap_from_df
from utils.run_r_gene_heatmap import generate_heatmap_from_df
from utils.run_r_gene_violin_plot import generate_violin_plot_from_df
from utils.run_r_gene_dot_plot import generate_dot_plot_from_df

# Load environment variables from .env file
load_dotenv()

# --- Data Caching ---
OPTIONS_CACHE = {} # Cache for dropdown options (clusters, subjects)
GENE_LIST_CACHE = {} # Cache for gene lists
JOBS_IN_PROGRESS = set()
REFRESHED_JOBS = set()
GENE_REFRESH_FLAGS = {}
s3_client = boto3.client('s3')

# --- REFRESH FLAG HELPERS ---
def set_refresh_flag(dataset_prefix, status):
    """Set refresh flag state: idle | running | ready."""
    GENE_REFRESH_FLAGS[dataset_prefix] = {
        "status": status,
        "timestamp": time.time(),
    }
    log_progress(f"🔖 Flag updated for {dataset_prefix}: {status}")

def get_refresh_flag(dataset_prefix):
    """Return current refresh state (idle | running | ready)."""
    return GENE_REFRESH_FLAGS.get(dataset_prefix, {}).get("status", "idle")

# --- GENE UNIVERSE LOADER ---
def get_available_gene_universe(dataset_prefix, bucket_name=None, force_s3=False):
    # 1. Load Bucket from Env if not provided
    if not bucket_name:
        bucket_name = os.getenv("S3_BUCKET_URI")
    
    if not bucket_name:
        print("Error: S3_BUCKET_URI not set in .env or passed as argument.")
        return []
    
    # 2. Construct paths
    s3_key = f"Joe/HSV_Dashboard_py/DataWarehouse/GEX/{dataset_prefix}_avail_genelist.json"
    local_path = f"DataWarehouse/GEX/{dataset_prefix}_avail_genelist.json"

    # 3. Check Local File (Only if force_s3 is False)
    if not force_s3 and os.path.exists(local_path):
        print(f"Loading from LOCAL: {local_path}")
        with open(local_path, "r") as f:
            return json.load(f)

    # 4. Load from S3
    try:
        # Clean the bucket name if it contains "s3://"
        if bucket_name.startswith("s3://"):
            bucket_name = urlparse(bucket_name).netloc

        print(f"Loading from S3: {bucket_name}/{s3_key}")
        
        obj = s3_client.get_object(Bucket=bucket_name, Key=s3_key)
        return json.loads(obj["Body"].read().decode('utf-8'))

    except Exception as e:
        print(f"Error loading from S3: {e}")
        return []
    
# --- DATA LOADING HELPERS --- 
def get_dataset_options(dataset_prefix, bucket_name=None, force_s3=False):
    """Loads UMAP data to get unique clusters and subjects."""
    # 0. Check Memory Cache First
    if dataset_prefix in OPTIONS_CACHE and not force_s3:
        return OPTIONS_CACHE[dataset_prefix]
    
    # 1. Load Bucket from Env if not provided
    if not bucket_name:
        bucket_name = os.getenv("S3_BUCKET_URI")
    
    if not bucket_name:
        print("Error: S3_BUCKET_URI not set")
        return {"clusters": [], "subjects": []}
    
    # 2. Construct paths
    s3_key = f"Joe/HSV_Dashboard_py/DataWarehouse/UMAP/{dataset_prefix}_umap_data.parquet"
    local_path = f"DataWarehouse/UMAP/{dataset_prefix}_umap_data.parquet"

    df = None

    try:
        # 3. Check local file
        if not force_s3 and os.path.exists(local_path):
            df = pd.read_parquet(local_path)
        # 4. Load from S3
        else:
        # Clean the bucket name if it contains "s3://"
            if bucket_name.startswith("s3://"):
                bucket_name = urlparse(bucket_name).netloc

            print(f"Loading from S3: {bucket_name}/{s3_key}")
            
            obj = s3_client.get_object(Bucket=bucket_name, Key=s3_key)
            df = pd.read_parquet(BytesIO(obj['Body'].read())) 

        # 5. Extract options and cache
        if df is not None:
            OPTIONS_CACHE[dataset_prefix] = {
                "clusters": sorted(df['CellType_Level3'].unique()),
                "subjects": sorted(df['Subject'].unique())
            }
            return OPTIONS_CACHE[dataset_prefix]
        
    except FileNotFoundError as e:
            print(f"Error loading dataset options for {dataset_prefix}: {e}")
            return {"clusters": [], "subjects": []}
    
    return OPTIONS_CACHE[dataset_prefix]

def get_gene_list(dataset_prefix, bucket_name=None, force_s3=False):
    """Loads the pre-computed list of available genes from the JSON index."""
    # 1. Check Memory Cache First
    if dataset_prefix in GENE_LIST_CACHE and not force_s3:
        return GENE_LIST_CACHE[dataset_prefix]
    
    # 2. Load Bucket from Env if not provided
    if not bucket_name:
        bucket_name = os.getenv("S3_BUCKET_URI")
    
    if not bucket_name:
        print("Error: S3_BUCKET_URI not set")
        return {"clusters": [], "subjects": []}
    
    # 3. Construct paths
    s3_key = f"Joe/HSV_Dashboard_py/DataWarehouse/GEX/{dataset_prefix}_gex_genes.json"
    local_path = f"DataWarehouse/GEX/{dataset_prefix}_gex_genes.json"

    data = []

    try:
        # 4. Check local file
        if not force_s3 and os.path.exists(local_path):
            with open(local_path, 'r') as f:
                data = json.load(f)
        # 5. Load from S3
        else:
            if bucket_name:
                if bucket_name.startswith("s3://"):
                    bucket_name = urlparse(bucket_name).netloc

                print(f"Loading from S3: {bucket_name}/{s3_key}")
                obj = s3_client.get_object(Bucket=bucket_name, Key=s3_key)
                data = json.loads(obj["Body"].read().decode('utf-8'))
            
        # 6. Cache and return
        if data:
            GENE_LIST_CACHE[dataset_prefix] = data
            return data
    except Exception as e:
        print(f"Error loading gene list for {dataset_prefix}: {e}")
        return []
    return []

def check_genes_availability(dataset_prefix, genes, bucket_name=None, force_s3=False):
    """
    Check which genes are available (Locally or S3) vs which need to be pulled.
    """
    # 1. Force a fresh check by clearing the cache for this dataset.
    # This ensures we pick up any recent updates from the R script or S3.
    if dataset_prefix in GENE_LIST_CACHE:
        del GENE_LIST_CACHE[dataset_prefix]

    # 2. Load the gene list using the robust helper we just created.
    # This handles the .env loading, S3 fallback, and path logic automatically.
    available_list = get_gene_list(dataset_prefix, bucket_name=bucket_name, force_s3=force_s3)
    available_genes = set(available_list)

    # 3. Handle case where list is completely missing (empty)
    if not available_genes:
        print(f"Gene index not found for {dataset_prefix}. Treating all as missing.")
        return [], genes

    # 4. Calculate Delta (Found vs Missing)
    found = [g for g in genes if g in available_genes]
    missing = [g for g in genes if g not in available_genes]

    return found, missing   

def build_ordered_gene_list(dropdown_genes, typed_genes):
    """Combine dropdown + typed → ordered unique uppercase genes."""
    genes = dropdown_genes or []
    if typed_genes:
        for g in typed_genes.replace(",", " ").split():
            g = g.strip().upper()
            if g and g not in genes:
                genes.append(g)
    return genes

# --- NON-BLOCKING BACKGROUND FETCH (NEW) ---
def log_progress(message):
    """Simple timestamped logger for backend progress tracking."""
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] {message}")

def run_precompute_r_async(dataset_prefix, genes_to_extract):
    """
    Fire-and-forget call to the R script to update parquet and JSON for missing genes.
    Adds detailed print logging to track backend progress.
    """
    if not genes_to_extract:
        return

    key = (dataset_prefix, tuple(sorted(genes_to_extract)))
    if key in JOBS_IN_PROGRESS:
        log_progress(f"Background R job already running for {dataset_prefix}: {genes_to_extract}")
        return

    JOBS_IN_PROGRESS.add(key)

    r_script_path = "utils/precompute_exp_run.R"
    env = os.environ.copy()
    env["EXTRACT_PREFIX"] = dataset_prefix
    env["EXTRACT_GENES"] = ",".join(genes_to_extract)

    log_progress(f"Launching async R process for dataset '{dataset_prefix}' to pull genes: {genes_to_extract}")
    try:
        process = subprocess.Popen(
            ["Rscript", r_script_path, "--vanilla"],
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        log_progress(f"PID {process.pid} started. Monitoring in background...")

        # background watcher to log progress
        import threading

        def monitor_process(proc, job_key):
            try:
                for line in iter(proc.stdout.readline, ""):
                    if line.strip():
                        log_progress(f"[R stdout] {line.strip()}")
                for err_line in iter(proc.stderr.readline, ""):
                    if err_line.strip():
                        log_progress(f"[R stderr] {err_line.strip()}")
                proc.wait()
                if proc.returncode == 0:
                    log_progress(f"R job completed successfully for {dataset_prefix}: {genes_to_extract}")

                    if dataset_prefix in GENE_LIST_CACHE:
                        del GENE_LIST_CACHE[dataset_prefix]
                    log_progress(f"🧹 Cleared gene list cache for {dataset_prefix}")

                    set_refresh_flag(dataset_prefix, "ready") 
                else:
                    log_progress(f"R job exited with code {proc.returncode} for {dataset_prefix}")
            except Exception as e:
                log_progress(f"Error monitoring R process: {e}")
            finally:
                JOBS_IN_PROGRESS.discard(job_key)
                log_progress(f"Job key {job_key} removed from in-progress list.")

        threading.Thread(target=monitor_process, args=(process, key), daemon=True).start()

    except Exception as e:
        log_progress(f"Failed to launch async R job: {e}")
        JOBS_IN_PROGRESS.discard(key)

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
                log_progress("🧊 No local genes available; skipping immediate plot generation.")
    
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
