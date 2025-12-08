import os
import boto3
import json
import tempfile
from urllib.parse import urlparse

# Global S3 client
s3_client = boto3.client("s3")

def get_bucket():
    bucket = os.getenv("S3_BUCKET_URI")
    if not bucket:
        raise ValueError("S3_BUCKET_URI not set")

    if bucket.startswith("s3://"):
        bucket = urlparse(bucket).netloc

    return bucket

def load_local_or_s3_parquet(local_path, bucket, key):
    # local available?
    if os.path.exists(local_path):
        print(f"USING LOCAL FILE: {local_path}")
        return local_path

    # fallback to S3
    print(f"USING S3 FILE: s3://{bucket}/{key}")
    obj = s3_client.get_object(Bucket=bucket, Key=key)
    data_bytes = obj["Body"].read()

    with tempfile.NamedTemporaryFile(suffix=".parquet", delete=False) as tmp:
        tmp.write(data_bytes)
        return tmp.name
    
def load_s3_umap(dataset_prefix, force_s3=True):
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    local_path = os.path.join(project_root, "DataWarehouse/UMAP", f"{dataset_prefix}_umap_data.parquet")

    # Local mode (allowed when force_s3=False)
    if not force_s3 and os.path.exists(local_path):
        print(f"USING LOCAL UMAP: {local_path}")
        return local_path

    # S3 mode
    bucket = get_bucket()
    key = f"Joe/HSV_Dashboard_py/DataWarehouse/UMAP/{dataset_prefix}_umap_data.parquet"
    print("DEBUG S3 KEY REQUEST:", key)

    obj = s3_client.get_object(Bucket=bucket, Key=key)
    data_bytes = obj["Body"].read()

    # Write to temp file → ALWAYS RETURN VALID FILE PATH
    with tempfile.NamedTemporaryFile(suffix=".parquet", delete=False) as tmp:
        tmp.write(data_bytes)
        return tmp.name

def load_s3_stats_cluster_status(dataset_prefix, force_s3=True):
    bucket = get_bucket()
    local_path = f"DataWarehouse/Stat/{dataset_prefix}_stats_cluster_status.parquet"
    key = f"Joe/HSV_Dashboard_py/DataWarehouse/Stat/{dataset_prefix}_stats_cluster_status.parquet"
    return load_local_or_s3_parquet(local_path, bucket, key)


def load_s3_stats_cluster_sample(dataset_prefix, force_s3=True):
    bucket = get_bucket()
    local_path = f"DataWarehouse/Stat/{dataset_prefix}_stats_cluster_sample.parquet"
    key = f"Joe/HSV_Dashboard_py/DataWarehouse/Stat/{dataset_prefix}_stats_cluster_sample.parquet"
    return load_local_or_s3_parquet(local_path, bucket, key)

def load_s3_stats_subject_status(dataset_prefix, force_s3=True):
    bucket = get_bucket()
    local_path = f"DataWarehouse/Stat/{dataset_prefix}_stats_subject_status.parquet"
    key = f"Joe/HSV_Dashboard_py/DataWarehouse/Stat/{dataset_prefix}_stats_subject_status.parquet"
    return load_local_or_s3_parquet(local_path, bucket, key)

def load_s3_colors(dataset_prefix, force_s3=True):
    bucket = get_bucket()
    local_path = f"DataWarehouse/Color/{dataset_prefix}_colors.json"
    key = f"Joe/HSV_Dashboard_py/DataWarehouse/Color/{dataset_prefix}_colors.json"

    # LOCAL file exists?
    if not force_s3 and os.path.exists(local_path):
        print(f"USING LOCAL COLOR JSON: {local_path}")
        text = open(local_path).read().strip()
        # Auto-fix Python dict
        if text.startswith("{'"):
            text = text.replace("'", '"')
        return json.loads(text)

    # Otherwise load from S3
    print("DEBUG S3 KEY REQUEST:", key)
    obj = s3_client.get_object(Bucket=bucket, Key=key)
    text = obj["Body"].read().decode("utf-8")
    # Auto-fix Python dict syntax if needed
    if text.startswith("{'"):
        text = text.replace("'", '"')

    return json.loads(text)
