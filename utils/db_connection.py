# utils/db_connection.py
import duckdb
import threading
import boto3
import os

_lock = threading.Lock()
_con = None

def get_duckdb():
    global _con
    with _lock:
        if _con is None or _con.closed:
            _con = duckdb.connect(database=':memory:')
        return _con

def configure_duckdb_s3(con):
    """
    Manually passes Boto3's active credentials to DuckDB.
    This ensures DuckDB works exactly like Boto3 (Local & ECS).
    """
    try:
        # 1. Install httpfs (Required for S3)
        con.execute("INSTALL httpfs; LOAD httpfs;")
        
        # 2. Fetch fresh credentials from Boto3
        # Boto3 handles looking at .env, ~/.aws/, or ECS Task Roles automatically.
        session = boto3.Session()
        creds = session.get_credentials()
        region = session.region_name or os.getenv("AWS_DEFAULT_REGION", "us-west-2")

        if not creds:
            print("⚠️ Warning: No AWS credentials found via Boto3.")
            return

        # Force Boto3 to refresh if creds are almost expired
        creds = creds.get_frozen_credentials()

        # 3. Pass keys to DuckDB
        con.execute(f"SET s3_region='{region}';")
        con.execute(f"SET s3_access_key_id='{creds.access_key}';")
        con.execute(f"SET s3_secret_access_key='{creds.secret_key}';")
        
        if creds.token:
            con.execute(f"SET s3_session_token='{creds.token}';")

    except Exception as e:
        print(f"❌ Failed to configure DuckDB S3: {e}")

def dict_to_r_vector(py_dict):
    """
    Converts Python dict: {'Cluster A': '#FF0000', 'Cluster B': '#00FF00'} 
    To R string: 'c("Cluster A"="#FF0000", "Cluster B"="#00FF00")'
    """
    if not py_dict or not isinstance(py_dict, dict):
        return "NULL"
    
    # Format as "key"="value", safely quoting both sides
    items = [f'"{str(k)}"="{str(v)}"' for k, v in py_dict.items()]
    return f"c({', '.join(items)})"
