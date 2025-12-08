# utils/run_r_spatial_deconvo.py

import os
import boto3
from urllib.parse import urlparse
from dotenv import load_dotenv
import plotly.graph_objects as go
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

# Load env variables
load_dotenv()

# Initialize s3 cilent
s3_client = boto3.client('s3')

# -----------------------------------------------------------------------------
# Load Data Once at Startup
# -----------------------------------------------------------------------------
# The path to your actual R data file
DATA_FILE = "DataWarehouse/Visium/Data__Visium_Deconvolution_Interactive.Rdata" 

if not os.path.exists(DATA_FILE):
    print(f"Local file missing: {DATA_FILE}. Attempting S3 download...")
    
    # 2. Attempt Download from S3
    try:
        # Ensure the local directory exists
        os.makedirs(os.path.dirname(DATA_FILE), exist_ok=True)
        
        # Get Bucket Config
        bucket_uri = os.getenv("S3_BUCKET_URI")
        if not bucket_uri:
            raise ValueError("S3_BUCKET_URI not set in .env")
            
        bucket_name = urlparse(bucket_uri).netloc
        # Construct Key
        s3_key = f"Joe/HSV_Dashboard_py/{DATA_FILE}"
        
        print(f"Downloading from s3://{bucket_name}/{s3_key}...")
        
        s3 = boto3.client('s3')
        s3.download_file(bucket_name, s3_key, DATA_FILE)
        print("Download successful.")
        
    except Exception as e:
        print(f"Failed to download S3 file: {e}")

# 3. Final Validation
if not os.path.exists(DATA_FILE):
    raise FileNotFoundError(f"CRITICAL: Data file could not be found locally or downloaded from S3: {DATA_FILE}")
        
# Normalize slashes for R
data_file_r = DATA_FILE.replace("\\", "/")

# Load the Rdata file into the rpy2 environment
ro.r(f'''
suppressMessages({{
  suppressWarnings({{
    load("{data_file_r}")
  }})
}})
''')

# Convert the R objects to Python objects
with localconverter(ro.default_converter + pandas2ri.converter):
    Res4Plot_A1 = ro.r['Res4Plot_A1']
    cell_type_names_A1 = list(ro.r['cell_type_names_A1'])
    # R named vectors need to be converted to Python dicts
    colors_r = ro.r['colors']
    colors = {name: colors_r.rx(name)[0] for name in colors_r.names}

# -----------------------------------------------------------------------------
def run_r_spatial_deconvo():
    """
    Generate Visium deconvolution pie plots from pre-loaded data.
    """
    fig = go.Figure()

    # Compute axis limits
    x_min = Res4Plot_A1["imagecol"].min() - 55
    x_max = Res4Plot_A1["imagecol"].max() + 55
    y_min = Res4Plot_A1["imagerow"].max() + 55
    y_max = Res4Plot_A1["imagerow"].min() - 55

    pie_size = 0.02  # relative size of each pie

    # Add pie charts for each row
    for idx, row in Res4Plot_A1.iterrows():
        pie_values = row[cell_type_names_A1].values
        pie_labels = cell_type_names_A1

        # Only include non-zero slices
        non_zero_idx = [i for i, v in enumerate(pie_values) if v > 0]
        if non_zero_idx:
            fig.add_trace(go.Pie(
                labels=[pie_labels[i] for i in non_zero_idx],
                values=[pie_values[i] for i in non_zero_idx],
                marker=dict(colors=[colors.get(pie_labels[i]) for i in non_zero_idx]),
                textinfo="none",
                hoverinfo="label+percent",
                domain=dict(
                    x=[(row.imagecol - x_min)/(x_max - x_min) - pie_size,
                       (row.imagecol - x_min)/(x_max - x_min) + pie_size],
                    y=[(row.imagerow - y_min)/(y_max - y_min) - pie_size,
                       (row.imagerow - y_min)/(y_max - y_min) + pie_size]
                ),
                showlegend=False
            ))

    # Layout
    fig.update_layout(
        title="",
        xaxis=dict(title="Image Column", range=[x_min, x_max], zeroline=False, showgrid=False),
        yaxis=dict(title="Image Row", range=[y_max, y_min], zeroline=False, showgrid=False, autorange="reversed"),
        showlegend=True,
        width=1000,
        height=800,
        margin=dict(t=5, b=20, l=20, r=20),
        plot_bgcolor='white',
        paper_bgcolor='white'
    )
    # Ensure axes are visible
    fig.update_xaxes(showline=True, linewidth=1, linecolor='black', mirror=True)
    fig.update_yaxes(showline=True, linewidth=1, linecolor='black', mirror=True)

    return fig
