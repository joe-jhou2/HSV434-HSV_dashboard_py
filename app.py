# app.py
from dash import Dash, html, dcc, Input, Output
import dash_bootstrap_components as dbc
from pages.intro_page import intro_layout
from pages.main_page import main_layout
from pages.tabs import scrnaseq_cluster_tab
from pages.tabs import scrnaseq_gene_tab
from pages.tabs import visium_spatial_tab
from pages.tabs import visium_deconv_tab

# ----------------------------------------
app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP],
           suppress_callback_exceptions=True)

server = app.server

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

# ----------------------------------------
if __name__ == "__main__":
    app.run(debug=True)