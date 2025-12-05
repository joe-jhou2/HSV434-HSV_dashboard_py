# utils/run_gex_data_loader.py
import os
import glob
import time
import duckdb
import json
import pandas as pd
import boto3
from urllib.parse import urlparse
from dotenv import load_dotenv
from utils.db_connection import configure_duckdb_s3

# Load env
load_dotenv()
s3_client = boto3.client('s3')

def load_filtered_gex_data(dataset_prefix, genes=None, clusters=None, subjects=None, bucket_name=None, force_s3=False):
    """
    Load AND JOIN filtered GEX data from multiple Parquet files
    using pure Python and DuckDB.
    """
    # Start timing
    start_time = time.time()

    # --- 1. Setup Paths & Mode (Local vs S3) ---
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

    # S3 Path Definitions
    bucket_uri = os.getenv("S3_BUCKET_URI")

    if bucket_uri.startswith("s3://"):
        actual_bucket = urlparse(bucket_uri).netloc

    s3_prefix = "Joe/HSV_Dashboard_py/DataWarehouse/GEX"

    # Local Path Definitions
    local_gex_dir = os.path.join(project_root, "DataWarehouse", "GEX")
    local_core_path = os.path.join(local_gex_dir, f"{dataset_prefix}_gex_core.parquet")

    # State Variables
    use_s3 = False
    core_path = ""
    ext_files = []
    color_map = []

    # Determine if we should use local files or S3
    if not force_s3 and os.path.exists(local_core_path):
        print(f"📂 Using LOCAL GEX files from: {local_gex_dir}")
        use_s3 = False
        core_path = local_core_path

        # Local: Find extension files using Glob
        all_files_pattern = os.path.join(local_gex_dir, f"{dataset_prefix}_gex_*.parquet")
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
             print("❌ Error: No local file and S3_BUCKET_URI is missing.")
             return pd.DataFrame(), {}

        print(f"☁️ Using S3 GEX files from bucket: {actual_bucket}")
        use_s3 = True

        # S3: Core file path
        core_path = f"s3://{actual_bucket}/{s3_prefix}/{dataset_prefix}_gex_core.parquet"

        try: 
            # S3: List objects to find core and extension files
            response = s3_client.list_objects_v2(Bucket=actual_bucket,
                                                 Prefix=f"{s3_prefix}/{dataset_prefix}_gex_*.parquet")
            all_files = [obj['Key'] 
                         for obj in response.get('Contents', []) 
                         if obj['Key'].startswith(f"{s3_prefix}/{dataset_prefix}_gex_*.parquet")]
            core_key = f"{s3_prefix}/{dataset_prefix}_gex_core.parquet"
            ext_files = [key for key in all_files if key != core_key]

        except Exception as e:
            print(f"❌ Error listing S3 files: {e}")
            return pd.DataFrame(), {}
            
        # S3: Load Colors
        color_key = f"Joe/HSV_Dashboard_py/DataWarehouse/Color/{dataset_prefix}_colors.json"
        try:
            obj = s3_client.get_object(Bucket=actual_bucket, Key=color_key)
            color_map = json.loads(obj['Body'].read().decode('utf-8'))
        except Exception as e:
            print(f"❌ Error loading color file from S3: {e}")

    # --- 2. Get Core Schema (to identify duplicate columns) ---
    # Connect to an in-memory DuckDB instance
    con = duckdb.connect()

    if use_s3:
        configure_duckdb_s3(con)

    try:
        core_schema_df = con.execute(f"DESCRIBE SELECT * FROM read_parquet('{core_path}')").df()
        core_cols = set(core_schema_df['column_name'])
        print(f"Core file has {len(core_cols)} columns (metadata, genes, etc.)")
    except Exception as e:
        print(f"❌ Error reading schema for core file {core_path}: {e}")
        con.close()
        return pd.DataFrame(), color_map
    
    # --- 3. Dynamically Build the SQL Query ---
    # Base of the query
    from_clause = f"FROM read_parquet('{core_path}') AS core"
    join_clauses = []
    
    # Keep track of all columns we've seen. Start with the core columns.
    all_seen_cols = set(core_cols)
    # Map a column name to its table alias (e.g., 'CD3E' -> 't1')
    col_to_table_map = {col: 'core' for col in core_cols}

    # Loop through extension files to build the JOIN clauses
    for i, file_path in enumerate(ext_files):
        alias = f't{i}'
        try:
            ext_schema_df = con.execute(f"DESCRIBE SELECT * FROM read_parquet('{file_path}')").df()
            ext_cols = set(ext_schema_df['column_name'])
            
            # This is your column pruning logic:
            new_cols = ext_cols - all_seen_cols
            
            if "Barcode" not in ext_cols:
                print(f"⚠️ Skipping {os.path.basename(file_path)}: No 'Barcode' column.")
                continue
                
            if not new_cols:
                print(f"ℹ️ Skipping {os.path.basename(file_path)}: No new columns found.")
                continue

            print(f"   └─ Joining {os.path.basename(file_path)} (alias {alias}) for {len(new_cols)} new columns.")

            # Add the new columns to our maps
            all_seen_cols.update(new_cols)
            for col in new_cols:
                col_to_table_map[col] = alias

            # Add the JOIN clause for this file
            join_clauses.append(
                f"LEFT JOIN read_parquet('{file_path}') AS {alias} ON core.Barcode = {alias}.Barcode"
            )

        except Exception as e:
            print(f"❌ Error reading schema for {file_path}: {e}")
            continue

    # --- 4. Build the final SELECT and WHERE clauses ---
    # Define the "metadata" columns we always want
    required_cols = {"Barcode", "UMAP_1", "UMAP_2", "CellType_Level3", "Subject", "Status"}
    
    # Add the specific genes the user requested
    cols_to_select = set(required_cols)
    if genes:
        gene_list = genes if isinstance(genes, list) else [genes]
        cols_to_select.update(gene_list)
    else:
        # If no genes specified, select ALL available columns
        cols_to_select.update(all_seen_cols)

    # Build the final SELECT list, using the correct table alias for each column
    final_select_list = []
    missing_cols = set()
    for col in cols_to_select:
        if col in col_to_table_map:
            table_alias = col_to_table_map[col]
            # Use dot notation: core."UMAP_1", t0."CD3E", etc.
            final_select_list.append(f'{table_alias}."{col}"')
        else:
            if col in (gene_list if genes else []): # Only warn for user-requested genes
                missing_cols.add(col)

    if missing_cols:
        print(f"⚠️ Warning: Requested genes not found in any file: {missing_cols}")

    # Build WHERE clause (filtering on the 'core' table's metadata)
    where_clauses = ["1=1"]
    if clusters:
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
        print(f"✅ Joined {len(ext_files) + 1} files into {df.shape[0]:,} rows × {df.shape[1]} cols")
        return df, color_map
    except Exception as e:
        print(f"❌ DuckDB Query Failed: {e}")
        return pd.DataFrame(), color_map
    finally:
        elapsed_time = time.time() - start_time
        print(f"⏱️ load_filtered_gex_data() completed in {elapsed_time:.2f} seconds.")
        con.close()
