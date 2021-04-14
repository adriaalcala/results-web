import dash
import dash_core_components as dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output
import dash_html_components as html
import dash_table


from flask import request, Flask

from dashboard.pinaweb_results.callback import update_graph, update_table, layout, get_info_data
from texts import GRAPH_DESCRIPTION, TABLE_DESCRIPTION
# ======== #
# DASH APP #
# ======== #
# Dash uses the Flask web framework
# For the Dash deployment, we need to access the Flask application instance
FONT_AWESOME = "https://use.fontawesome.com/releases/v5.10.2/css/all.css"


app = dash.Dash(
    __name__,
    url_base_pathname='/util-aligner-results/',
    external_stylesheets=[dbc.themes.BOOTSTRAP, FONT_AWESOME],
)

app.layout =html.Div(
    [
        dcc.Location(id='url', refresh=False), html.H1('PINAWEB RESULTS', id='H1-text', style={'textAlign': 'center'}),
        dcc.Tabs([
            dcc.Tab(label='Comparison', children=[
                dcc.Markdown(GRAPH_DESCRIPTION, style={'margin': '60px'}),
                html.Div(id='page-content', children=[  
                    html.Div(id='concensus-plot'),
                html.P(TABLE_DESCRIPTION, style={'margin': '60px'}),
                    html.Div([
                        html.P('Page size: ', id='page-size-text', style={'margin-left': '60px'}),  
                        dcc.Input(id='table-concensus-plot-page-count', type='number', min=1, max=100000, value=20, style={'margin': '6px'})
                    ],
                    style={'display': 'flex'}
                    ),
                    html.Div(id='div-dash-table', children=dash_table.DataTable(
                            id='table-concensus-plot',
                            columns=[],
                            page_current=0,
                            page_size=20,
                            page_action='custom',
                            filter_action='custom',
                            filter_query='',
                            sort_action='custom',
                            sort_mode='multi',
                            sort_by=[],
                            style_table={'overflowY': 'auto', 'overflowX': 'auto', 'width': '80%', 'margin_left': '10px'},
                            style_cell={
                                'height': 'auto',
                                # all three widths are needed
                                'minWidth': '10px', 'width': '10px', 'maxWidth': '10px',
                                'whiteSpace': 'normal'
                            }

                        )
                    ),
                ])
            ]),
            dcc.Tab(label='Alignments info', id="info"),

        ])
    ]
)


@app.callback(
    [Output('table-concensus-plot', "data"),
    Output('table-concensus-plot', "columns"),
    Output('concensus-plot', "children"),
    Output('info', "children")],
    [Input('table-concensus-plot', "page_current"),
    Input('table-concensus-plot', "page_size"),
    Input('table-concensus-plot', "sort_by"),
    Input('table-concensus-plot', "filter_query")])
def update_consensus_table(page_current, page_size, sort_by, filter_query):
    job_id = request.headers['Referer'].split('/')[-1]
    rows, columns = update_table(job_id, page_current, page_size, sort_by, filter_query)
    info_data = get_info_data(job_id)
    return rows, columns, update_graph(rows), info_data

@app.callback(
    Output('table-concensus-plot', "page_size"),
     [Input('table-concensus-plot-page-count', 'value'),])
def update_page_count(page_count_value):
    if page_count_value is None:
        return 20
    return page_count_value


# ======= #
# RUN APP #
# ======= #
server = app.server
app.config['suppress_callback_exceptions'] = True

# Run Dash on a Public IP
if __name__ == "__main__":
    app.run_server(host="0.0.0.0", port=int("5000"), debug=False, dev_tools_hot_reload = False, use_reloader=True)
