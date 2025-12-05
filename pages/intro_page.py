# pages/intro_page.py

from dash import html, Input, Output, callback
import dash_bootstrap_components as dbc

def intro_layout():
    return dbc.Container([
        html.Link(rel="stylesheet", href="styles.css"),

        html.Div(className="FHCRC-image",
                 children=html.Img(src="/assets/images/FHRC_logo.png")),

        html.Div([
            html.H1("Welcome to HSV Genital Skin Biopsy scRNAseq Explorer",
                    style={"textAlign": "center"}),
            html.H2("Introduction", style={"textAlign": "center"}),
            html.P(
                "This study explores the composition of immune cells in genital skin tissue from HSV-infected patients using scRNA-seq data.",
                style={"fontSize": "20px", "textAlign": "center"}
            ),
            html.P(
                "Our application focuses on examining the dynamics of cell types (clusters) across different disease statuses and visualizing gene expression patterns within each cluster and disease status.",
                style={"fontSize": "20px", "textAlign": "center"}
            ),
            html.Br(),
            html.P("Click below to proceed to the Explorer",
                   style={"fontSize": "18px", "color": "blue", "textAlign": "center"}),
            html.Div(style={"textAlign": "center"}, children=[
                dbc.Button("Proceed to Explorer",
                           id="goToMain",
                           color="primary",
                           style={"fontSize": "18px", "padding": "10px 20px"})
            ]),
            html.Br(),
            html.Div(html.Img(src="/assets/images/HSV.png"),
                     className="Virus-image"),
            html.Div(html.Img(src="/assets/images/Work_Flow_ver2.png",
                              style={"width": "120%", "height": "auto"}),
                     className="WorkFlow-image")
        ])
    ], fluid=True)

# --- Redirect to /main when button clicked ---
@callback(
        Output("url", "pathname"),
        Input("goToMain", "n_clicks"),
        prevent_initial_call=True
        )
def go_to_main(n_clicks):
    if n_clicks:
        print("Redirecting to /main")
        return "/main"

