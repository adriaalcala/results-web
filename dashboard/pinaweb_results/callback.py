import plotly.figure_factory as ff
import dash_core_components as dcc
import dash_html_components as html
import dash_table
import pandas as pd
import requests
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
from collections import Counter

df = pd.DataFrame() #pd.read_csv('dashboard/pinaweb_results/data/sample_2.tsv', delimiter='\t')


def layout(job_id="5f609561ceadad3aecd73bd5"):
    global df
    print(job_id)
    aligner = None
    try:
        file_id = requests.get(f'https://biocom.uib.es/util-aligner/v2/comparison/{job_id}').json()['files']['joined_tsv']
    except:
        data = requests.get(f'https://biocom.uib.es/util-aligner/v2/alignment/{job_id}').json()
        file_id = data['files']['alignment_tsv']
        aligner = data['aligner']
    df = pd.read_csv(f'https://biocom.uib.es/util-aligner/v2/file/{file_id}', delimiter='\t')
    if aligner:
        df = df.rename(columns={df.columns[1]: f"{df.columns[1]}_{aligner}"})
    species_1 = df.columns[0].split('_')[-1]
    species_2 = df.columns[1].split('_')[-2]
    url = 'http://ltimbigd2.uib.es:11080/db/stringdb/items/species/select'
    req_headers = {'Accept': 'text/tab-separated-values'}

    response = requests.post(url=url, json={'columns': ['species_id', 'official_name'], 'filter':{'species_id':[int(species_1), int(species_2)]}}, headers=req_headers)
    species_dict = dict([i.split('\t') for i in response.text.split('\n')[:-1]])
    
    url = 'http://ltimbigd2.uib.es:11080/db/stringdb/items/proteins/select'
    response = requests.post(url=url, json={'columns':['protein_external_id','preferred_name'], 'filter':{'species_id':[int(species_1), int(species_2)]}}, headers=req_headers)
    proteins_dict = dict([i.split('\t') for i in response.text.split('\n')[:-1]])

    df = df.rename(columns={df.columns[0]: species_dict[species_1]})
    df[f"{species_dict[species_1]}\npreferred name"] = df.apply(lambda row: proteins_dict.get(row[0]), axis=1)
    num_columns = df.shape[1]
    for i in range(1, num_columns-1):
        aligner_name = df.columns[i].split('_')[-1]
        df = df.rename(columns={df.columns[i]: f"{species_dict[species_2]}_{aligner_name}"})
        df[f"{species_dict[species_2]}_{aligner_name}\npreferred name"] = df.apply(lambda row: proteins_dict.get(row[i]), axis=1)
    if not aligner:
        df['Concensus'] = df.apply(lambda row: max(Counter(row).values())/(df.shape[1]-1), axis=1)
    
    return [html.Div(id='concensus-plot'),
            'Page size: ',
            dcc.Input(
            id='table-concensus-plot-page-count',
            type='number',
            min=1,
            max=100000,
            value=200
            ),
            html.Div(dash_table.DataTable(
                id='table-concensus-plot',
                columns=[{"name": i, "id": i} for i in sorted(df.columns)],
                page_current=0,
                page_size=200,
                page_action='custom',
                filter_action='custom',
                filter_query='',
                sort_action='custom',
                sort_mode='multi',
                sort_by=[],
                style_table={'overflowY': 'auto', 'overflowX': 'auto'},
                style_cell={
                    'height': 'auto',
                    # all three widths are needed
                    'minWidth': '10px', 'width': '10px', 'maxWidth': '10px',
                    'whiteSpace': 'normal'
                }
            ),
            ),
    ]

def update_graph(rows):
    all_target_proteins = set()
    alignment_keys = [k for k in rows[0].keys() if '_' in k and '\n' not in k]
    origin_key = [k for k in rows[0].keys() if k not in alignment_keys and '\n' not in k][0]
    for row in rows:
        all_target_proteins = all_target_proteins.union({row[k] for k in alignment_keys})
    all_target_proteins = list(all_target_proteins)
    origin = []
    aligners = [i.split('_')[-1] for i in alignment_keys]
    print(aligners)
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
        hover.append([ 'Origin: ' + origin[idx] + '<br>' 'Target: ' + i + '<br>' + 'Aligner: ' + aligners[x]
                        for idx, i in enumerate(alignments[x])])
    if len(aligners) == 1:
        return [] 
    fig = ff.create_annotated_heatmap(
        z=alignments_num, annotation_text=alignments_text, text=hover, hoverinfo='text', y=aligners, x=origin,
        colorscale=[[0, "black"], [0.5, "red"], [1.0, "green"]], showscale=True, zmin=-1, zmax=1
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


def update_table(page_current, page_size, sort_by, filter):
    filtering_expressions = filter.split(' && ')
    dff = df
    for filter_part in filtering_expressions:
        col_name, operator, filter_value = split_filter_part(filter_part)

        if operator in ('eq', 'ne', 'lt', 'le', 'gt', 'ge'):
            # these operators match pandas series operator method names
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

    return dff.iloc[ 
        page_current*page_size: (page_current + 1)*page_size
    ].to_dict('records')
