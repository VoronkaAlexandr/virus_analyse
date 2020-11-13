import json
import os
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import requests
import csv
import yaml
from collections import Counter
import numpy as np

root_path = Path('..')
data_path = root_path / "data"
raw_data_path = data_path /'raw'
processed_data_path = data_path / 'processed'
figures = root_path / 'figures'

def process_node(node, parent):
    # node['branch_attrs', 'name', 'node_attrs']
    name = node['name']
    # num_descendants = descendants_count(node)
    attrs = {}
    mutations = {}
    p = -1
    children_count = bfs_count_children(node)
    for node_attr in node['node_attrs']:
        val = node['node_attrs'][node_attr]
        if isinstance(val, dict) and 'value' in val:
            # print(node_attr, val['value'])
            attrs[node_attr] = val['value']

    if 'div' in node['node_attrs']:
        attrs['div'] = node['node_attrs']['div']

    if 'branch_attrs' in node and 'mutations' in node['branch_attrs']:
        mutations = node['branch_attrs']

    if parent is not None:
        p = node['node_attrs']['tree_id']

    return (node['node_attrs']['tree_id'], name, p, children_count, attrs, mutations)

def bfs_tree(tree):
    q = []
    q.append(tree)
    idx = 0
    while len(q) > 0:
        n = q.pop(0)
        n['node_attrs']['tree_id'] = idx
        idx += 1
        if 'children' in n:
            for child in n['children']:
                q.append(child)

def get_children(tree, nodes, parent):
    nodes.append(process_node(tree, parent))
    if 'children' in tree:
        for child in tree['children']:
            get_children(child, nodes, tree)
            # nodes.append(process_node(child, tree))

def bfs_count_children(tree):
    cnt = 0
    q = []
    q.append(tree)

    while len(q) > 0:
        n = q.pop(0)
        cnt += 1
        if 'children' in n:
            for child in n['children']:
                q.append(child)
    return cnt - 1

def extract_data_from_json(fn):
    with open(os.path.join('..', 'data', 'raw', '{}.json'.format(fn))) as f:
        data = json.load(f)
    tree = data['tree']

    bfs_tree(tree)
    nodes_list = []

    get_children(tree, nodes_list, None)
    k1 = list({'age',
               'author',
               'clade_membership',
               'country',
               'division',
               'gisaid_epi_isl',
               'host',
               'location',
               'num_date',
               'originating_lab',
               'recency',
               'region',
               'submitting_lab'}) + ['div']  # ,'labels', 'mutations']
    with open('{}.tsv'.format(fn), 'w') as f:
        f.write('\t'.join(['idx', 'parent', 'name', 'child_count'] + k1 + ['labels', 'mutations']) + os.linesep)
        for node_item in nodes_list:
            node_id = node_item[0]
            name = node_item[1]
            parent = node_item[2]
            child_count = node_item[3]
            _attrs = node_item[4]
            mut = node_item[-1]
            mut_part = ['', '']
            if 'labels' in mut:
                mut_part[0] = str(mut['labels'])
            if 'mutations' in mut and len(mut['mutations']) > 0:
                mut_part[1] = str(mut['mutations'])

            attr_line = [_attrs.get(key, '') for key in k1]

            #         line = '\t'.join([node_item[0]]
            #                          + [str(node_item[1][key]) for key in k1 ]
            #                          +mut_part
            #                         )

            line = [node_id, parent, name, child_count] + attr_line + mut_part
            f.write('\t'.join([str(itm) for itm in line]) + os.linesep)
    df = pd.read_csv('{}.tsv'.format(fn), sep='\t')
    df.rename(columns={'div': 'divergence'}, inplace=True)

    result_csv_dir = processed_data_path / fn
    if not result_csv_dir.exists():
        result_csv_dir.mkdir()

    df.to_csv(result_csv_dir / '{}.tsv'.format(fn), sep='\t', index=None, index_label=None)
    plt.figure(figsize=(10, 8))
    sns.distplot(df['divergence'])
    plt.title('Divergence distribution')
    plt.grid()
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.savefig(result_csv_dir / 'divergence_distribution_{}.svg'.format(fn))
    plt.figure(figsize=(10, 8))
    plt.scatter(df['divergence'], df['num_date'] - df['num_date'].min())
    plt.xlabel('Divergence')
    plt.ylabel('Time')
    plt.title('Time-Divegence')
    plt.grid()
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.savefig(result_csv_dir / 'joint_divergence_time_{}.svg'.format(fn))

def download_json(name, fn):
    url =  name
    resp = requests.get(url)

    if resp.status_code == 200:
        with open(raw_data_path / fn , 'w') as f:
            f.write(resp.text)
        return True

    return False

def mut_spec():
    with open(processed_data_path / 'covid41020/covid41020.tsv') as tsvfile:
        tsv_table = pd.read_csv(tsvfile, sep='\t')
        from_to = []
        where_mut = []
        for line in tsv_table['mutations'].dropna():
            my_mutations = yaml.load(line)['nuc']
            for single_mut in my_mutations:
                if ((single_mut[-1]) != '-') and ((single_mut[0]) != '-'):
                    from_to.append((single_mut[0], single_mut[-1]))
                    where_mut.append(single_mut[1:-1])
        print(Counter(from_to))
        print(where_mut)
        covid_mutations = pd.DataFrame(
            {'From': list(map(lambda x: x[0], Counter(from_to).keys())),
             'To': list(map(lambda x: x[-1], Counter(from_to).keys())),
             'probability': Counter(from_to).values()})
        covid_mutations.to_csv(processed_data_path / 'covid41020/covid41020_mut_spectr.tsv', sep='\t')
        file = open(processed_data_path / "covid41020/where_was_mutations.txt", "w")
        file.write(str(where_mut))
        file.close()


if __name__ == '__main__':
     base_url = 'https://nextstrain.org/charon/getDataset?prefix=/ncov/global'
     if download_json(base_url, "covid41020.json"):
         if extract_data_from_json('covid41020'):
             mut_spec()
