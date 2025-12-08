import os
import glob
import duckdb
import datetime

# --- Configuration ---
DATA_PATH = "DataWarehouse/Pert"
TARGETS = ["tcell", "myeloid"] # Prefixes to process
LOG_FILE = "DataWarehouse/logs/pert_merge.log"

def log(msg):
    """Logs a message to stdout and to the log file."""
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line)
    os.makedirs(os.path.dirname(LOG_FILE), exist_ok=True)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")

def safe_path(p):
    """Helper to ensure file paths have forward slashes for DuckDB SQL."""
    return p.replace(os.sep, '/')

JOIN_KEYS = ["Subject", "CellType_Level3", "Status"]

def merge_pert_files_duckdb(prefix):
    """
    Merge all parquet files for a dataset into a single core parquet
    using DuckDB for memory-efficient out-of-core merging.
    """
    log(f"--- Checking {prefix} with DuckDB ---")
    
    # --- 1. Identify Core File ---
    core_path = os.path.join(DATA_PATH, f"{prefix}_pert_core.parquet")
    if not os.path.exists(core_path):
        log(f"No core file found for {prefix}. Skipping.")
        return

    # --- 2. Find all *other* parquet files ---
    pattern = os.path.join(DATA_PATH, f"{prefix}_pert_*.parquet")
    all_files = glob.glob(pattern)
    ext_files = [
        f for f in all_files if not os.path.normpath(f).endswith(os.path.normpath(core_path))
    ]

    if not ext_files:
        log(f"No new parquet files found for {prefix}.")
        return

    # Connect to an in-memory DuckDB
    con = duckdb.connect()
    
    # Create a "view" or temp table for the core file
    con.execute(f"CREATE VIEW core AS SELECT * FROM read_parquet('{safe_path(core_path)}')")
    
    try:
        # Get core columns to identify what's "new"
        core_schema_df = con.execute("DESCRIBE core").df()
        core_cols = set(core_schema_df['column_name'])
        log(f"Core file loaded ({len(core_cols)} cols). Checking {len(ext_files)} extension files...")

        join_clauses_sql = []
        all_new_cols_sql = []
        files_merged_successfully = []

        # --- 3. Build a single, massive SQL query ---
        for i, f_path in enumerate(ext_files):
            alias = f't{i}'
            try:
                ext_schema_df = con.execute(f"DESCRIBE SELECT * FROM read_parquet('{safe_path(f_path)}')").df()
                ext_cols = set(ext_schema_df['column_name'])

                if not all(k in ext_cols for k in JOIN_KEYS):
                    log(f"Skipping {os.path.basename(f_path)} — missing join keys.")
                    continue
                
                # Find only columns that are NOT in the core file
                new_cols = [
                    c for c in ext_cols
                    if c not in core_cols and c not in JOIN_KEYS
                ]

                if not new_cols:
                    log(f"{os.path.basename(f_path)}: no new columns.")
                    files_merged_successfully.append(f_path)
                    continue


                # Add the JOIN clause for this file
                join_condition = " AND ".join([f"core.{k} = {alias}.{k}" for k in JOIN_KEYS])
                
                join_clauses_sql.append(
                    f"LEFT JOIN read_parquet('{safe_path(f_path)}') AS {alias} ON {join_condition}"
                )
                
                # Add the new columns to the final SELECT statement
                all_new_cols_sql.extend([f'{alias}."{c}"' for c in new_cols])
                
                # If we get here, the file is good to merge and delete
                files_merged_successfully.append(f_path)
                log(f"Adding {len(new_cols)} new columns from {os.path.basename(f_path)}")

            except Exception as e:
                log(f"Error inspecting {os.path.basename(f_path)}: {e}")
        
        if not all_new_cols_sql:
            log(f"No new columns found in any files for {prefix}.")
            # We still run cleanup for any files that had no new columns
            if files_merged_successfully:
                log("--- Cleaning up processed files ---")
                for f in files_merged_successfully:
                    try:
                        os.remove(f)
                        log(f"Deleted {os.path.basename(f)}")
                    except Exception as e:
                        log(f"Failed to delete {os.path.basename(f)}: {e}")
            return

        # --- 4. Execute the merge to a temporary file ---
        temp_output_path = core_path + ".tmp"
        
        final_select_sql = "core.*, " + ", ".join(all_new_cols_sql)
        
        query = f"""
        COPY (
            SELECT {final_select_sql}
            FROM core
            {' '.join(join_clauses_sql)}
        ) TO '{safe_path(temp_output_path)}' (FORMAT 'parquet', CODEC 'zstd');
        """
        
        log(f"Executing merge query for {len(files_merged_successfully)} files...")
        con.execute(query)
        con.close() # Release the connection

        # --- 5. Atomic Replace and Cleanup ---
        # If SQL succeeds, atomically replace the old core file with the new one
        os.replace(temp_output_path, core_path)
        log(f"Successfully saved new core file: {os.path.basename(core_path)}")

        log("--- Cleaning up merged files ---")
        for f in files_merged_successfully:
            try:
                os.remove(f)
                log(f"Deleted {os.path.basename(f)}")
            except Exception as e:
                log(f"Failed to delete {os.path.basename(f)}: {e}")

    except Exception as e:
        log(f"CRITICAL MERGE FAILED for {prefix}: {e}")
        con.close()
        # Clean up the temp file if it exists
        if os.path.exists(temp_output_path):
            os.remove(temp_output_path)
            log("Cleaned up temp file. Core file is untouched.")

def daily_merge():
    log("--- Daily PERT Merge Job STARTING---")
    for prefix in TARGETS:
        try:
            merge_pert_files_duckdb(prefix)
        except Exception as e:
            log(f"Unhandled error for {prefix}: {e}")
    log("---Daily PERT Merge Job FINISHED---")

# --- For testing, run it once immediately ---
log("Running one-time merge for testing...")
daily_merge()
log("Test run complete.")