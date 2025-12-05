# app.py
from dash import Dash, html, dcc, Input, Output
import dash_bootstrap_components as dbc
from pages.intro_page import intro_layout
from pages.main_page import main_layout
from pages.tabs import scrnaseq_cluster_tab
from pages.tabs import scrnaseq_gene_tab
from pages.tabs import visium_spatial_tab
from pages.tabs import visium_deconv_tab

# from apscheduler.schedulers.background import BackgroundScheduler
# from utils.auto_gex_parquet_merger import daily_merge

# ----------------------------------------
app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP],
           suppress_callback_exceptions=True)

app.layout = html.Div([
    dcc.Location(id="url", refresh=False),
    html.Div(id="page-content")
])

# Register callbacks and pass the thread-safe function
scrnaseq_cluster_tab.register_callbacks(app)
scrnaseq_gene_tab.register_callbacks(app)
visium_spatial_tab.register_callbacks(app)
visium_deconv_tab.register_callbacks(app)

# -------- Page Routing --------
@app.callback(Output("page-content", "children"),
              Input("url", "pathname"))
def display_page(pathname):
    if pathname in [None, "/"]:
        return intro_layout()
    elif pathname == "/main":
        return main_layout()
    else:
        return html.H3("404: Page not found")

# --- Logging setup ---
# LOG_DIR = "DataWarehouse/logs"
# os.makedirs(LOG_DIR, exist_ok=True)
# LOG_FILE = os.path.join(LOG_DIR, "scheduler.log")

# def scheduler_log(msg):
#     ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
#     line = f"[{ts}] {msg}"
#     print(line)
#     with open(LOG_FILE, "a") as f:
#         f.write(line + "\n")

# ----------------------------------------
if __name__ == "__main__":
    app.run(debug=True)
    # Start scheduler AFTER Dash server
    # scheduler = BackgroundScheduler()
    # scheduler.add_job(daily_merge, 'cron', hour=0, minute=0, next_run_time=datetime.now() + timedelta(days=1))
    # scheduler.start()
    # scheduler_log("🕒 Scheduler started — daily_merge runs every day at 00:00 (midnight).")
