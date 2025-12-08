import os
import time
import json
import boto3
import datetime
import subprocess
from urllib.parse import urlparse

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
