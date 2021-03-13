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
#df = pd.DataFrame() #pd.read_csv('dashboard/pinaweb_results/data/sample_2.tsv', delimiter='\t')


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
    species_1 = df.columns[0].split('_')[-1]
    species_2 = df.columns[1].split('_')[-2]
    url = 'http://ltimbigd2.uib.es:11080/db/stringdb/items/species/select'
    req_headers = {'Accept': 'text/tab-separated-values'}

    response = requests.post(url=url, json={'columns': ['species_id', 'official_name'], 'filter':{'species_id':[int(species_1), int(species_2)]}}, headers=req_headers)
    species_dict = dict([i.split('\t') for i in response.text.split('\n')[:-1]])
    
    url = 'http://ltimbigd2.uib.es:11080/db/stringdb/items/proteins/select'
    response = requests.post(url=url, json={'columns':['protein_external_id','preferred_name', 'annotation'], 'filter':{'species_id':[int(species_1), int(species_2)]}}, headers=req_headers)
    proteins_dict = dict([i.split('\t')[:2] for i in response.text.split('\n')[:-1]])

    df = df.rename(columns={df.columns[0]: species_dict[species_1]})
    df[f"{species_dict[species_1]} preferred name"] = df.apply(lambda row: proteins_dict.get(row[0]), axis=1)
    columns_to_drop = [f"{species_dict[species_1]}"]
    num_columns = df.shape[1]
    for i in range(1, num_columns-1):
        aligner_name = df.columns[i].split('_')[-1]
        df = df.rename(columns={df.columns[i]: f"{species_dict[species_2]}_{aligner_name}"})
        df[f"{species_dict[species_2]}_{aligner_name} preferred name"] = df.apply(lambda row: proteins_dict.get(row[i]), axis=1)
        columns_to_drop.append(f"{species_dict[species_2]}_{aligner_name}")
    
    df = df.drop(columns_to_drop, axis=1)
    if not aligner:
        df['Concensus'] = df.apply(lambda row: float("{:.2f}".format(max(Counter(row).values())/(df.shape[1]-1))), axis=1)
    df_dicts[job_id] = df
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
    aligners = [i.split('_')[-1].split()[0] for i in alignment_keys]
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
        colorscale=[[0, "grey"], [0.5, "red"], [1.0, "green"]], showscale=False, zmin=-1, zmax=1
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
    columns = [{"name": i.replace('preferred name', '').replace('_', ''), "id": i} for i in sorted(dff.columns)]
    return dff2, columns
