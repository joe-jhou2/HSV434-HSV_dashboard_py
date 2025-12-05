# pages/tabs/visium_spatial_tab.py

import dash_bootstrap_components as dbc
from dash import html, dcc, Input, Output, State
from utils.run_r_spatial_gene import run_r_spatial_gene

# -----------------------------------------------------------------------------
# Layout for Visium Spatial Tab
# -----------------------------------------------------------------------------
spatial_tab_layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            html.Label("Enter gene name: (separate multiple genes with commas)"),
            dbc.Input(id="gene_name_input", type="text", placeholder="e.g., CD3D, CD8A, GZMB")
        ], width=4),
        dbc.Col([
            dbc.Button("Run Analysis", id="run_gene_spatial", color="primary")
        ], width=2)
    ]),
    html.Br(),
    dcc.Loading(
        id="spatial_loading",
        type="circle",
        children=html.Div(id="dynamic_spatial_plot_ui")
    ),
    # Optional modal for loading feedback
    dbc.Modal(
        [
            dbc.ModalHeader(dbc.ModalTitle("Analyzing Data")),
            dbc.ModalBody("Please wait while the spatial plot is generated...")
        ],
        id="spatial_loading_modal",
        is_open=False,
        backdrop="static",
        keyboard=False
    )
], fluid=True)

# -----------------------------------------------------------------------------
# Callback
# -----------------------------------------------------------------------------
def register_callbacks(app):
    @app.callback(
        Output("dynamic_spatial_plot_ui", "children"),
        Output("spatial_loading_modal", "is_open"),
        Input("run_gene_spatial", "n_clicks"),
        State("gene_name_input", "value"),
        prevent_initial_call=True
    )
    def update_spatial_plot(n_clicks, gene_name):
        if not gene_name:
            # Return empty content and keep modal closed
            return html.Div("Please enter a gene name and click 'Run Analysis'."), False

        plot_component, modal_state = run_r_spatial_gene(gene_name)
        return plot_component, modal_state
