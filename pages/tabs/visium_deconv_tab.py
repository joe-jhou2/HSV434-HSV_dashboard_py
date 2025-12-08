# pages/tabs/visium_deconv_tab.py

from dash import html, dcc, Input, Output
import dash_bootstrap_components as dbc
from utils.run_r_spatial_deconvo import run_r_spatial_deconvo

# -----------------------------------------------------------------------------
# Layout for Visium Deconvolution Tab
# -----------------------------------------------------------------------------
deconv_tab_layout = dbc.Container([
    # html.H3("Cell Type Deconvolution", style={"textAlign": "center"}),
    html.Br(),
    dbc.Row([
        dbc.Col([
            dbc.Button("Run Analysis", id="run_deconvo", color="primary"),
        ], width=2),
        dbc.Col([
            dcc.Loading(
                id="deconv_loading",
                type="circle",
                children=html.Div(id="dynamic_pie_plot_ui")
            )
        ], width=10)
    ])
], fluid=True)

# -----------------------------------------------------------------------------
# Callbacks
# -----------------------------------------------------------------------------
def register_callbacks(app):
    @app.callback(
        Output("dynamic_pie_plot_ui", "children"),
        Input("run_deconvo", "n_clicks"),
        prevent_initial_call=True
    )
    def update_deconv_plot(n_clicks):
        if not n_clicks:
            return html.Div("Click 'Run Analysis' to generate plots.", style={"color": "gray"})

        try:
            # Call the Python function that generates the figure
            fig = run_r_spatial_deconvo()  # gene_name optional for future use
            return dcc.Graph(figure=fig)
        except Exception as e:
            return html.Pre(f"Error generating deconvolution plot:\n{e}")
