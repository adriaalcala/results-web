import dash
import dash_core_components as dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output
import dash_html_components as html

from flask import request, Flask

from dashboard.pinaweb_results.callback import update_graph, update_table, layout
# ======== #
# DASH APP #
# ======== #
# Dash uses the Flask web framework
# For the Dash deployment, we need to access the Flask application instance
FONT_AWESOME = "https://use.fontawesome.com/releases/v5.10.2/css/all.css"

server = Flask(__name__)

app = dash.Dash(
    __name__,
    routes_pathname_prefix="/",
    external_stylesheets=[dbc.themes.BOOTSTRAP, FONT_AWESOME],
    server=server
)

app.layout =html.Div([dcc.Location(id='url', refresh=False), html.H1('PINAWEB RESULTS', style={'text-align': 'center'}), html.Div(id='page-content')])


@app.callback(
    Output('page-content', 'children'),
    [Input('url', 'pathname')])
def display_page(pathname):
    print(pathname)
    return layout(job_id=pathname.split('/')[-1])


@app.callback(
    Output('concensus-plot', "children"),
    [Input('table-concensus-plot', "data")])
def plot_concensus(rows):
    return update_graph(rows)


@app.callback(
    Output('table-concensus-plot', "data"),
    [Input('table-concensus-plot', "page_current"),
    Input('table-concensus-plot', "page_size"),
    Input('table-concensus-plot', "sort_by"),
    Input('table-concensus-plot', "filter_query")])
def update_consensus_table(page_current, page_size, sort_by, filter_query):
    print(page_current, page_size)
    return update_table(page_current, page_size, sort_by, filter_query)

@app.callback(
    Output('table-concensus-plot', "page_size"),
     [Input('table-concensus-plot-page-count', 'value'),])
def update_page_count(page_count_value):
    if page_count_value is None:
        return 20
    return page_count_value

@app.callback([Output("protein-images", "children")],
              [Input("dropdown-proteines", "value")],
              )
def get_protein_images(protein):
    return None


# ======= #
# RUN APP #
# ======= #
server = app.server
app.config['suppress_callback_exceptions'] = True

# Run Dash on a Public IP
if __name__ == "__main__":
    app.run_server(host="0.0.0.0", port=int("8080"), debug=False, use_reloader=True)