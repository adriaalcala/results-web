import requests
import pandas as pd
from collections import Counter


jobs = [
    '6075b48df9ef164d857914a6', '6075b5d6f9ef164d857914bf', '6075b620f9ef164d857914d6', '6075b6f3f9ef164d857914ed',
    '6075b725f9ef164d85791504', '6075b767f9ef164d8579151b', '6075b81ff9ef164d85791534', '6075bb65f9ef164d8579154b', 
    '6075bc5bf9ef164d85791562', '6075be1ff9ef164d85791579'
]
nets = {
    10090: 'Mus musculus',
    6239: 'Caenorhabditis elegans',
    4932: 'Saccharomyces cerevisiae',
    9606: 'Homo sapiens',
    7227: 'Drosophila melanogaster'
}


def compute_consensus(row):
    n = len(row)
    clean_row = [i for i in row if type(i) == str]
    return max(Counter(clean_row).values())/n

results = []
for job_id in jobs:
    dd = requests.get(f'https://biocom.uib.es/util-aligner/v2/comparison/{job_id}').json()
    file_id = dd['files']['joined_tsv']

    df = pd.read_csv(f'https://biocom.uib.es/util-aligner/v2/file/{file_id}', delimiter='\t')

    df['Consensus'] = df.apply(compute_consensus, axis=1)

    results.append([nets[dd['net1']['species_id']], nets[dd['net2']['species_id']], sum(df['Consensus'])/len(df['Consensus'])])

print(results)
print(max(results, key=lambda x: x[2]))
print(min(results, key=lambda x: x[2]))