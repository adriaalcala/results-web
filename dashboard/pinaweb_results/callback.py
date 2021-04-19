import plotly.figure_factory as ff
import dash_core_components as dcc
import dash_html_components as html
import dash_table
import pandas as pd
import requests
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
from collections import Counter

df_dicts = {} 
response_dict = {}
#df = pd.DataFrame() #pd.read_csv('dashboard/pinaweb_results/data/sample_2.tsv', delimiter='\t')

ALIGNERS = {
    'alignet': 'AligNet',
    'hubalign': 'HubAlign',
    'l-graal': 'L-GRAAL', 
    'pinalog': 'PINALOG',
    'spinal': 'SPINAL'
}


def compute_consensus(row):
    n = len(row) - 1
    clean_row = [i for i in row[1:] if type(i) == str]
    if not clean_row:
        return 0
    return max(Counter(clean_row).values())/n


def layout(job_id="5f609561ceadad3aecd73bd5"):
    global df_dicts
    if job_id in df_dicts:
        return df_dicts[job_id]
    aligner = None
    try:
        print('file_url', f'https://biocom.uib.es/util-aligner/v2/comparison/{job_id}')
        file_id = requests.get(f'https://biocom.uib.es/util-aligner/v2/comparison/{job_id}').json()['files']['joined_tsv']
        print('file_id', file_id)
    except:
        data = requests.get(f'https://biocom.uib.es/util-aligner/v2/alignment/{job_id}').json()
        file_id = data['files']['alignment_tsv']
        aligner = data['aligner']
    df = pd.read_csv(f'https://biocom.uib.es/util-aligner/v2/file/{file_id}', delimiter='\t')
    if aligner:
        df = df.rename(columns={df.columns[1]: f"{df.columns[1]}_{aligner}"})
    if 'custom' in df.columns[0]:
        species_1 = '-1'
        proteins_species_1 = df.columns[0].split('_')[-1]
    else:
        species_1 = df.columns[0].split('_')[-1]
        proteins_species_1 = df.columns[0].split('_')[-1]
    if 'custom' in df.columns[1]:
        species_2 = '-1'
        proteins_species_2 = df.columns[1].split('_')[-2]
    else:
        species_2 = df.columns[1].split('_')[-2]
        proteins_species_2 = df.columns[1].split('_')[-2]

    url = 'http://ltimbigd2.uib.es:11080/db/stringdb/items/species/select'
    req_headers = {'Accept': 'text/tab-separated-values'}
    response = requests.post(url=url, json={'columns': ['species_id', 'official_name'], 'filter':{'species_id':[int(species_1), int(species_2)]}}, headers=req_headers)
    species_dict = dict([i.split('\t') for i in response.text.split('\n')[:-1]])
    if species_1 not in species_dict:
        species_dict[species_1] = 'Custom'
    if species_2 not in species_dict:
        species_dict[species_2] = 'Custom'
    url = 'http://ltimbigd2.uib.es:11080/db/stringdb/items/proteins/select'
    response = requests.post(url=url, json={'columns':['protein_external_id','preferred_name', 'annotation'], 'filter':{'species_id':[int(proteins_species_1), int(proteins_species_2)]}}, headers=req_headers)
    proteins_dict = dict([i.split('\t')[:2] for i in response.text.split('\n')[:-1]])

    df = df.rename(columns={df.columns[0]: species_dict[species_1]})
    df[f"{species_dict[species_1]} preferred name"] = df.apply(lambda row: proteins_dict.get(row[0]), axis=1)
    columns_to_drop = [f"{species_dict[species_1]}"]
    num_columns = df.shape[1]
    print('df-1', df)
    for i in range(1, num_columns-1):
        aligner_name_raw  = df.columns[i].split('_')[-1]
        aligner_name = ALIGNERS.get(aligner_name_raw, aligner_name_raw)
        df = df.rename(columns={df.columns[i]: f"{species_dict[species_2]}_{aligner_name}"})
        df[f"{species_dict[species_2]}_{aligner_name} preferred name"] = df.apply(lambda row: proteins_dict.get(row[i]), axis=1)
        print()
        columns_to_drop.append(f"{species_dict[species_2]}_{aligner_name}")
    
    df = df.drop(columns_to_drop, axis=1)
    if not aligner:
        df['Consensus'] = df.apply(compute_consensus, axis=1)
    if 'Consensus' in df:
        df = df[df['Consensus'] > 0]
    df_dicts[job_id] = df
    print('df-2', df)
    columns = [{"name": i.replace('preferred name', ''), "id": i} for i in sorted(df.columns) if 'preferred' in i or 'Concensus' in i]
    return df


def update_graph(rows):
    all_target_proteins = set()
    if not rows:
        return dcc.Graph()
    alignment_keys = [k for k in rows[0].keys() if '_' in k]
    origin_key = [k for k in rows[0].keys() if k not in alignment_keys][0]
    for row in rows:
        all_target_proteins = all_target_proteins.union({row[k] for k in alignment_keys})
    all_target_proteins = list(all_target_proteins)
    origin = []
    aligners_raw = [i.split('_')[-1].split()[0] for i in alignment_keys]
    aligners = [ALIGNERS.get(i, i) for i in aligners_raw]
    alignments_text = [['']*len(rows) for _ in range(len(alignment_keys))]
    alignments = [[] for _ in range(len(alignment_keys))]
    alignments_num = [[] for _ in range(len(alignment_keys))]
    for row in rows:
        origin.append(row[origin_key])
        targets = Counter(row.values())
        for i in range(len(aligners)):
            if row[alignment_keys[i]]:
                #alignments_num[i].append(all_target_proteins.index(row[alignment_keys[i]]))
                alignments_num[i].append(targets[row[alignment_keys[i]]]/len(aligners))
                alignments[i].append(row[alignment_keys[i]])
            else:
                alignments_num[i].append(-1)
                alignments[i].append('-')

    hover=[]
    for x in range(len(aligners)):
        hover.append([ 'Origin: ' + str(origin[idx]) + '<br>' 'Target: ' + str(i) + '<br>' + 'Aligner: ' + str(aligners[x])
                        for idx, i in enumerate(alignments[x])])
    if len(aligners) == 1:
        return [] 
    fig = ff.create_annotated_heatmap(
        z=alignments_num, annotation_text=alignments_text, text=hover, hoverinfo='text', y=aligners, x=origin,
        colorscale=[[0, "grey"], [0.5, "white"], [1.0, "green"]], showscale=False, zmin=-1, zmax=1
        )
    fig.update_xaxes(showticklabels=False)

    return dcc.Graph(figure=fig)

operators = [['ge ', '>='],
             ['le ', '<='],
             ['lt ', '<'],
             ['gt ', '>'],
             ['ne ', '!='],
             ['eq ', '='],
             ['contains '],
             ['datestartswith ']]


def split_filter_part(filter_part):
    for operator_type in operators:
        for operator in operator_type:
            if operator in filter_part:
                name_part, value_part = filter_part.split(operator, 1)
                name = name_part[name_part.find('{') + 1: name_part.rfind('}')]

                value_part = value_part.strip()
                v0 = value_part[0]
                if (v0 == value_part[-1] and v0 in ("'", '"', '`')):
                    value = value_part[1: -1].replace('\\' + v0, v0)
                else:
                    try:
                        value = float(value_part)
                    except ValueError:
                        value = value_part

                # word operators need spaces after them in the filter string,
                # but we don't want these later
                return name, operator_type[0].strip(), value

    return [None] * 3


def update_table(job_id, page_current, page_size, sort_by, col_filter):
    filtering_expressions = col_filter.split(' && ')
    
    dff = layout(job_id)
    for filter_part in filtering_expressions:
        col_name, operator, filter_value = split_filter_part(filter_part)
        if operator in ('eq', 'ne', 'lt', 'le', 'gt', 'ge'):
            # these operators match pandas series operator method names
            try:
                filter_value = float(filter_value)
            except:
                pass
            dff = dff.loc[getattr(dff[col_name], operator)(filter_value)]
        elif operator == 'contains':
            dff = dff.loc[dff[col_name].str.contains(filter_value, na=False)]
        elif operator == 'datestartswith':
            # this is a simplification of the front-end filtering logic,
            # only works with complete fields in standard format
            dff = dff.loc[dff[col_name].str.startswith(filter_value)]
    if len(sort_by):
        dff = dff.sort_values(
            [col['column_id'] for col in sort_by],
            ascending=[
                col['direction'] == 'asc'
                for col in sort_by
            ],
            inplace=False
        )
    dff2 = dff.iloc[ 
        page_current*page_size: (page_current + 1)*page_size
    ].to_dict('records')
    columns = [{"name": i.replace('preferred name', '').replace('_', ' '), "id": i} for i in dff.columns]
    return dff2, columns


def get_info_data(job_id):
    info = requests.get(f'https://biocom.uib.es/util-aligner/v2/comparison/{job_id}').json()
    if not info:
        info  = requests.get(f'https://biocom.uib.es/util-aligner/v2/alignment/{job_id}').json()
    species_1 = str(info['net1']['species_id'])
    species_2 = str(info['net2']['species_id'])
    url = 'http://ltimbigd2.uib.es:11080/db/stringdb/items/species/select'
    req_headers = {'Accept': 'text/tab-separated-values'}
    scores = {}
    run_time = {}
    alignments = {}
    suplementary_data = {}
    if 'aligners' not in info:
        info['aligners'] = [info,]
    for alignment in info.get('results_object_ids', [job_id]):
        print(alignment)
        data = requests.get(f'https://biocom.uib.es/util-aligner/v2/alignment/{alignment}').json()
        print(data.keys())
        if 'scores' in data:
            scores[data['aligner']] = {
                'ec': data['scores']['ec_data']['ec_score'],
                'fc_score_jaccard': data['scores']['fc_data']['fc_score_jaccard'],
                'fc_score_hrss_bma': data['scores']['fc_data']['fc_score_hrss_bma']
            }
        run_time[data['aligner']] = data['results']['run_time']
        print('files', data['files'])
        alignments[data['aligner']] = data['files'].get('alignment_tsv')
        suplementary_data[data['aligner']] = alignment
    response = requests.post(url=url, json={'columns': ['species_id', 'official_name'], 'filter':{'species_id':[int(species_1), int(species_2)]}}, headers=req_headers)
    species_dict = dict([i.split('\t') for i in response.text.split('\n')[:-1]])
    edge_types_1 = ', '.join(f"{score}:{value}" for score, value in info['net1']['score_thresholds'].items())
    edge_types_2 = ', '.join(f"{score}:{value}" for score, value in info['net2']['score_thresholds'].items())
    style = {'border': '1px solid #aaaaaa', 'text-align': 'center', 'padding' : '8px'}
    style_2 = {'border': '1px solid #aaaaaa', 'text-align': 'center', 'padding' : '8px', 'background-color': "#dddddd"}
    style_3 = {'border': '1px solid #aaaaaa', 'text-align': 'center', 'padding' : '8px', 'background-color': "#eeeeee"}
    aligners_header_row = [
                html.Tr(children=[
                    html.Th("Aligner", style=style),
                    html.Th("Alignment", style=style),
                    html.Th("Run Time", style=style),
                    html.Th("Suplementary Data", style=style)
                ], style=style)
    ]
    aligners_row = [html.Tr(
        children=[
            html.Th(ALIGNERS.get(a['aligner'], a['aligner']), style=style_3),
            html.Th(html.A("alignment", href=f"https://biocom.uib.es/util-aligner/v2/file/{alignments[a['aligner']]}"), style=style_2 if idx%2 else style),
            html.Th(f"""{int(run_time[a['aligner']])}\"""", style=style_2 if idx%2 else style),
            html.Th(html.A("Suplementary data", href=f"https://biocom.uib.es/util-aligner/v2/alignment/{suplementary_data[a['aligner']]}"), style=style_2 if idx%2 else style),
        ])
        for idx, a in enumerate(sorted(info['aligners'], key=lambda x: x['aligner']))
    ]
    if 'joined' in info['results']:
        comparison_row = [
            html.Tr(
                children=[
                    html.Th('Comparison', style=style_3),
                    html.Th(html.A("joined", href=f"https://biocom.uib.es/util-aligner/v2/file/{info['results']['joined']['file']}"), style=style_2 if len(info['aligners'])%2 else style),
                    html.Th('', style=style_2 if len(info['aligners'])%2 else style),
                    html.Th(html.A(f"Suplementary data", href=f"https://biocom.uib.es/util-aligner/v2/comparison/{job_id}"), style=style_2 if len(info['aligners'])%2 else style)
                ]
            )
        ]
    else:
        comparison_row = []
    response = [
        html.H2('Results overview', style={ 'padding' : '16px', 'text-align': 'left' }),
        html.Table(
            children=[
                html.Tr(children=[
                    html.Th("Job results identifier", style={'text-align': 'center', 'padding' : '8px'}),
                    html.Th(job_id, style={'text-align': 'center', 'padding' : '8px'})
                ], style={'text-align': 'center', 'padding' : '8px'})
            ],
            style={ 'padding' : '8px', 'margin': '10px' }
        ),
        html.Div( children=[
        html.H4('Aligners', style={'padding': '8px' }),
        html.Table(
            children = aligners_header_row + aligners_row + comparison_row, style={'border': '1px solid #aaaaaa', 'text-align': 'center', 'padding' : '8px', 'width': '45%'}
        ),
        ],style={'align-self': 'center', 'padding': '8px'}),
        html.Div(children=[
            html.Div(children=[
                    html.H4("Input Network", style={ 'padding' : '8px' }),
                    html.Table(
                        children=[
                            html.Tr(children=[
                                html.Th("Species:", style=style_2),
                                html.Th(species_dict.get(species_1, 'Custom'), style=style_2)
                            ], style=style_2),                
                            html.Tr(children=[
                                html.Th("Edge types", style=style),
                                html.Th(edge_types_1, style=style)
                            ], style=style),
                            html.Tr(children=[
                                html.Th("Number of Proteins:", style=style_2),
                                html.Th(info['net1_details']['n_vert'], style=style_2)
                            ], style=style_2),
                            html.Tr(children=[
                                    html.Th("Number of Interactions:", style=style),
                                    html.Th(info['net1_details']['n_edges'], style=style)
                            ], style=style)
                        ],
                        style={ 'padding' : '8px', 'width': '100%' }
                    ),
                ], style={ 'padding' : '8px', 'width': '45%' }
            ),
            html.Div(children=[
                html.H4("Output Network", style={ 'padding' : '8px' }),
                html.Table(
                    children=[
                        html.Tr(children=[
                            html.Th("Species:", style=style_2),
                            html.Th(species_dict.get(species_2, 'Custom'), style=style_2)
                        ], style=style_2),                
                        html.Tr(children=[
                            html.Th("Edge types", style=style),
                            html.Th(edge_types_2, style=style)
                        ], style=style),
                        html.Tr(children=[
                            html.Th("Number of Proteins:", style=style_2),
                            html.Th(info['net2_details']['n_vert'], style=style_2)
                        ], style=style_2),
                        html.Tr(children=[
                                html.Th("Number of Interactions:", style=style),
                                html.Th(info['net2_details']['n_edges'], style=style)
                        ], style=style)
                    ],
                    style={ 'padding' : '8px', 'width': '100%' }
                ),
            ], style={ 'padding' : '8px', 'width': '45%' })
            ],
            style={'display': 'flex', 'align-self': 'center'}
        ),
        html.Div(children=[
            dcc.Graph(
                id='ec-score',
                figure={
                    'data': [
                        {'x': [ALIGNERS.get(a['aligner'], a['aligner']) for a in info['aligners']], 'y': [scores.get(i['aligner'], {'ec': 0})['ec'] for i in info['aligners']], 'type': 'bar', 'name': 'EC Score'},
                    ],
                    'layout': {
                        'title': 'EC Score'
                    }
                },
                style={'width': '45%'}
            ),
            dcc.Graph(
                id='fc-score',
                figure={
                    'data': [
                        {'x':  [ALIGNERS.get(a['aligner'], a['aligner']) for a in info['aligners']], 'y': [scores.get(i['aligner'], {}).get('fc_score_jaccard', 0) for i in info['aligners']], 'type': 'bar', 'name': 'Jaccard Score'},
                        {'x':  [ALIGNERS.get(a['aligner'], a['aligner']) for a in info['aligners']], 'y': [scores.get(i['aligner'], {}).get('fc_score_hrss_bma', 0) for i in info['aligners']], 'type': 'bar', 'name': 'HRSS Score'},
                    ],
                    'layout': {
                        'title': 'FC Score'
                    }
                },
                style={'width': '45%'}
            ),
        ], style={'display':'flex', 'align-self': 'center', 'align-items': 'center', 'width': '100%'}
        )

    ]
    return response

