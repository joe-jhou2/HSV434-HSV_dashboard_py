# pages/main_page.py

from dash import html, dcc
import dash_bootstrap_components as dbc
from pages.tabs.scrnaseq_cluster_tab import cluster_tab_layout
from pages.tabs.scrnaseq_gene_tab import gene_tab_layout
from pages.tabs.visium_spatial_tab import spatial_tab_layout
from pages.tabs.visium_deconv_tab import deconv_tab_layout

def main_layout():
    return dbc.Container([
        dcc.Tabs(id="main-tabs", value="scrna", children=[
            dcc.Tab(label="scRNAseq", value="scrna", className="scRNAseq-tab-label", children=[
                html.Br(),
                dbc.Row([
                    dbc.Col([
                        html.H5("Select Dataset"),
                        dcc.Dropdown(
                            id="dataset_option",
                            options=[
                                {"label": "T cell dataset", "value": "tcell"},
                                {"label": "Myeloid cell dataset", "value": "myeloid"}
                            ],
                            placeholder="Select dataset",
                            value="myeloid"
                        ),
                        html.P("Select a dataset and run analysis in the tabs below.", className="text-muted small mt-2")
                    ], width=6),
                ]),
                html.Hr(),
                dcc.Tabs([
                    dcc.Tab(label="Cluster Discovery", 
                            value="cluster_tab", 
                            className="Cluster-tab-label",
                            children=cluster_tab_layout),
                    dcc.Tab(label="Gene Discovery", 
                            value="gene_tab", 
                            className="Gene-tab-label",
                            children=gene_tab_layout)
                ])
            ]),
            dcc.Tab(label="Visium", value="visium", className="Visium-tab-label", children=[
                html.Div(className="visium-image",
                 children=html.Img(src="/assets/images/visium.png")),
                html.Br(),
                dcc.Tabs([
                    dcc.Tab(label="Gene Spatial Discovery",
                            value="spatial_tab",
                            className="Spatial-gene-tab-label",
                            children=spatial_tab_layout),
                    dcc.Tab(label="Cell Type Deconvolution",
                            value="deconv_tab",
                            className="Deconv-tab-label",
                            children=deconv_tab_layout)
                ])
            ])
        ])
    ], fluid=True)
