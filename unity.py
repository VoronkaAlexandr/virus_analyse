import pandas as pd
from pathlib import Path

root_path = Path('..')
data_path = root_path / "data"
raw_data_path = data_path /'raw'
processed_data_path = data_path / 'processed'
figures = root_path / 'figures'
unity = pd.DataFrame()
with open('../data/processed/covid41020/virus_data.txt', 'r') as virus_data:
    viruses = virus_data.readlines()
    for line in viruses:
        virus_name = line.split()[0]
        df = pd.read_csv(processed_data_path / 'mut_spectr/{}.tsv'.format(virus_name), sep = '\t')
        unity = unity.append(df)
unity.to_csv(processed_data_path / 'mut_spectr/unity_all.tsv', sep = '\t')