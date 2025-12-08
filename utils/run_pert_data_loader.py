# utils/run_pert_data_loader.py

import os
import glob
import duckdb
import pandas as pd
import time
import json
import boto3
from urllib.parse import urlparse
from dotenv import load_dotenv
from utils.db_connection import configure_duckdb_s3

# Load env
load_dotenv()

# Initialize S3 client
s3_client = boto3.client('s3')

# --- Key Columns for this dataset ---
KEY_COLS = {"Subject", "CellType_Level3", "Status"}

def safe_path(p):
    """Helper to ensure file paths have forward slashes for DuckDB SQL."""
    return p.replace(os.sep, '/')

def read_schema_duckdb(path):
    con = duckdb.connect(database=':memory:')
    df = con.execute(f"DESCRIBE SELECT * FROM read_parquet('{path}')").df()
    con.close()
    return set(df["column_name"])

def load_filtered_pert_data(dataset_prefix, genes=None, clusters=None, subjects=None, bucket_name=None, force_s3=False):
    """
    Loads and JOINS filtered Pert data from multiple Parquet files
    using DuckDB and a 'core' file as the base.
    """
    # Start timing
    start_time = time.time()

    # --- 1. Find Core and Extension Files ---
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

    # S3 Path Definitions
    bucket_uri = os.getenv("S3_BUCKET_URI")

    if bucket_uri.startswith("s3://"):
        actual_bucket = urlparse(bucket_uri).netloc

    s3_prefix = "Joe/HSV_Dashboard_py/DataWarehouse/Pert"

    # Local Path Definitions
    local_pert_dir = os.path.join(project_root, "DataWarehouse", "Pert")
    local_core_path = os.path.join(local_pert_dir, f"{dataset_prefix}_gex_core.parquet")

    # State Variables
    use_s3 = False
    core_path = ""
    ext_files = []
    color_map = []

    # Determine if we should use local files or S3
    if not force_s3 and os.path.exists(local_core_path):
        print(f"Using LOCAL Pert files from: {local_pert_dir}")
        use_s3 = False
        core_path = local_core_path

        # Local: Find extension files using Glob
        all_files_pattern = os.path.join(local_pert_dir, f"{dataset_prefix}_pert_*.parquet")
        all_files = glob.glob(all_files_pattern)
        ext_files = [f for f in all_files if os.path.normpath(f) != os.path.normpath(local_core_path)]
        
        # Local: Load Colors
        color_path = os.path.join(project_root, "DataWarehouse/Color", f"{dataset_prefix}_colors.json")
        if os.path.exists(color_path):
            with open(color_path, 'r') as f:
                color_map = json.load(f)
    else:
        # Use S3
        if not actual_bucket:
             print("Error: No local file and S3_BUCKET_URI is missing.")
             return pd.DataFrame(), {}

        print(f"Using S3 Pert files from bucket: {actual_bucket}")
        use_s3 = True

        # S3: Core file path
        core_path = f"s3://{actual_bucket}/{s3_prefix}/{dataset_prefix}_pert_core.parquet"

        try: 
            # S3: List objects to find core and extension files
            response = s3_client.list_objects_v2(Bucket=actual_bucket, 
                                                 Prefix=f"{s3_prefix}/{dataset_prefix}_pert_*.parquet")
            all_files = [obj['Key'] 
                         for obj in response.get('Contents', []) 
                         if obj['Key'].startswith(f"{s3_prefix}/{dataset_prefix}_pert_*.parquet")]
            core_key = f"{s3_prefix}/{dataset_prefix}_pert_core.parquet"
            ext_files = [key for key in all_files if key != core_key]

        except Exception as e:
            print(f"Error listing S3 files: {e}")
            return pd.DataFrame(), {}
            
        # S3: Load Colors
        color_key = f"Joe/HSV_Dashboard_py/DataWarehouse/Color/{dataset_prefix}_colors.json"
        try:
            obj = s3_client.get_object(Bucket=actual_bucket, Key=color_key)
            color_map = json.loads(obj['Body'].read().decode('utf-8'))
        except Exception as e:
            print(f"Error loading color file from S3: {e}")

    # --- 2. Get Core Schema (to identify duplicate columns) ---
    # Connect to an in-memory DuckDB instance
    con = duckdb.connect()

    if use_s3:
        configure_duckdb_s3(con)
    
    try:
        core_schema_df = con.execute(f"DESCRIBE SELECT * FROM read_parquet('{safe_path(core_path)}')").df()
        core_cols = set(core_schema_df['column_name'])
        print(f"Core file has {len(core_cols)} columns.")
        
        # Check if core file has the required keys
        if not KEY_COLS.issubset(core_cols):
            raise ValueError(f"Core file {os.path.basename(core_path)} is missing one or more keys: {KEY_COLS - core_cols}")

    except Exception as e:
        print(f"Error reading schema for core file {core_path}: {e}")
        con.close()
        return pd.DataFrame(), color_map

    # --- 3. Dynamically Build the SQL Query ---
    # Base of the query
    from_clause = f"FROM read_parquet('{safe_path(core_path)}') AS core"
    join_clauses = []
    
    # Map a column name to its table alias (e.g., 'GENE_X' -> 't1')
    col_to_table_map = {col: 'core' for col in core_cols}
    gene_list = [] # Initialize empty list for later check

    if genes:
        gene_list = genes if isinstance(genes, list) else [genes]

    # Loop through extension files to build the JOIN clauses
    for i, file_path in enumerate(ext_files):
        alias = f't{i}'
        try:
            ext_schema_df = con.execute(f"DESCRIBE SELECT * FROM read_parquet('{safe_path(file_path)}')").df()
            ext_cols = set(ext_schema_df['column_name'])
            
            # Check for keys
            if not KEY_COLS.issubset(ext_cols):
                print(f"Skipping {os.path.basename(file_path)}: Missing one or more keys.")
                continue
            
            # Find all potential new gene columns
            potential_new_cols = ext_cols - core_cols - KEY_COLS
            
            # Figure out which of these we actually need to join
            if genes:
                # If user specified genes, only join if this file has one of them
                cols_to_join = potential_new_cols.intersection(gene_list)
                if not cols_to_join:
                    print(f"Skipping {os.path.basename(file_path)}: No requested genes found.")
                    continue
            else:
                # If user did NOT specify genes, join all new columns
                cols_to_join = potential_new_cols
                if not cols_to_join:
                    print(f"Skipping {os.path.basename(file_path)}: No new columns found.")
                    continue
            
            print(f"Joining {os.path.basename(file_path)} (alias {alias}) for {len(cols_to_join)} columns.")

            # Update maps for the columns we're joining
            for col in cols_to_join:
                col_to_table_map[col] = alias

            # Add the JOIN clause for this file
            join_clauses.append(
                f"LEFT JOIN read_parquet('{safe_path(file_path)}') AS {alias} ON "
                f"core.Subject = {alias}.Subject AND "
                f"core.CellType_Level3 = {alias}.CellType_Level3 AND "
                f"core.Status = {alias}.Status"
            )

        except Exception as e:
            print(f"Error processing file {os.path.basename(file_path)} (alias {alias}): {e}")
            continue

    # --- 4. Build the final SELECT and WHERE clauses ---
    # Define the "metadata" columns
    required_cols = {"Subject", "CellType_Level3", "Status"}
    
    cols_to_select = set(required_cols)
    if genes:
        cols_to_select.update(gene_list)
    else:
        # If no genes specified, select ALL available columns
        cols_to_select.update(col_to_table_map.keys())

    # Build the final SELECT list, using the correct table alias
    final_select_list = []
    missing_cols = set()
    
    for col in cols_to_select:
        if col in col_to_table_map:
            table_alias = col_to_table_map[col]
            final_select_list.append(f'{table_alias}."{col}"')
        else:
            if col in gene_list: # Only warn for user-requested genes
                missing_cols.add(col)

    if missing_cols:
        print(f"Warning: Requested genes not found in any file: {missing_cols}")

    # Build WHERE clause (filtering on the 'core' table's metadata)
    where_clauses = ["1=1"]
    if clusters: # 'clusters' maps to 'CellType_Level3'
        cluster_sql_list = ", ".join([f"'{c}'" for c in clusters])
        where_clauses.append(f'core."CellType_Level3" IN ({cluster_sql_list})')
    if subjects:
        subject_sql_list = ", ".join([f"'{s}'" for s in subjects])
        where_clauses.append(f'core."Subject" IN ({subject_sql_list})')

    # --- 5. Assemble and Execute Final Query ---
    final_sql = f"""
    SELECT
        {', '.join(final_select_list)}
    {from_clause}
    {' '.join(join_clauses)}
    WHERE
        {' AND '.join(where_clauses)}
    """

    print("\n--- SQL TO EXECUTE ---")
    print(final_sql)
    print("----------------------")
    
    try:
        df = con.execute(final_sql).df()
        print(f"Joined {len(join_clauses) + 1} files into {df.shape[0]:,} rows × {df.shape[1]} cols")
        return df, color_map
    except Exception as e:
        print(f"DuckDB Query Failed: {e}")
        return pd.DataFrame(), color_map
    finally:
        elapsed_time = time.time() - start_time
        print(f"load_filtered_pert_data() completed in {elapsed_time:.2f} seconds.")
        con.close()
