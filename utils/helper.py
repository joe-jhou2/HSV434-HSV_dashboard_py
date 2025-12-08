import os
import boto3
from io import BytesIO
import pandas as pd
from urllib.parse import urlparse

# --- Data Caching ---
s3_client = boto3.client('s3')
OPTIONS_CACHE = {} # Cache for dropdown options (clusters, subjects)

def get_dataset_options(dataset_prefix, bucket_name=None, force_s3=False):
    """Loads a dataset into the cache and returns its unique clusters and subjects."""
    # 0. Check Cache First (Fast return)
    if dataset_prefix in OPTIONS_CACHE and not force_s3:
        return OPTIONS_CACHE[dataset_prefix]
    
    # 1. Load Bucket from Env if not provided
    if not bucket_name:
        bucket_name = os.getenv("S3_BUCKET_URI")
    
    if not bucket_name:
        print("Error: S3_BUCKET_URI not set in .env or passed as argument.")
        return {"clusters": [], "subjects": []}
    
    # 2. Construct paths
    s3_key = f"Joe/HSV_Dashboard_py/DataWarehouse/UMAP/{dataset_prefix}_umap_data.parquet"
    local_path = f"DataWarehouse/UMAP/{dataset_prefix}_umap_data.parquet"

    df = None

    try:
        # 3. Logic Branch: Local vs S3
        if not force_s3 and os.path.exists(local_path):
            df = pd.read_parquet(local_path)
        # 4. Load from S3
        else:
            # Clean the bucket name
            if bucket_name.startswith("s3://"):
                bucket_name = urlparse(bucket_name).netloc

            print(f"Loading UMAP options from S3: {bucket_name}/{s3_key}")
                
            obj = s3_client.get_object(Bucket=bucket_name, Key=s3_key)
            df = pd.read_parquet(BytesIO(obj['Body'].read()))

        # 5. Process and Cache
        if df is not None:
            OPTIONS_CACHE[dataset_prefix] = {
                "clusters": sorted(df['CellType_Level3'].unique()),
                "subjects": sorted(df['Subject'].unique())
            }
            return OPTIONS_CACHE[dataset_prefix]

    except Exception as e:
            print(f"Error loading dataset options for {dataset_prefix}: {e}")
            return {"clusters": [], "subjects": []}

    return OPTIONS_CACHE[dataset_prefix]
