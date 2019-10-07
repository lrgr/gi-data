import numpy as np
import pandas as pd
import sys
import click

from os.path import join
from tqdm import tqdm

from GIData import GIData
from GIData.utils import sparsity


'''
The file contains the complete SGA genetic interaction dataset in a tab-delimited format with 
columns: [
    'Query Strain ID', 
    'Query allele name', 
    'Array Strain ID',
    'Array allele name', 
    'Arraytype/Temp', # Array Type (DMA or TSA) and Temperature (26°C or 30°C)
    'Genetic interaction score (ε)',
    'P-value', 
    'Query single mutant fitness (SMF)', 
    'Array SMF',
    'Double mutant fitness', 
    'Double mutant fitness standard deviation']

Array Type (DMA or TSA) and Temperature (26°C or 30°C). Indicates the kind of array used:
    DMA{T} := non essential deletion mutant array
    TSA{T}:= essential TS allele array
where T is either 26 or 30

Possible Array Types:
DMA26
DMA30
TSA26
TSA30
'''

def id_to_name(s):
    return s.split('_')[0]

def load_gis(fp):
    # Renamed columns with shorter names...
    columns = [
        'query_strain_ID', 
        'query_allele_name', 
        'array_strain_ID',
        'array_allele_name', 
        'array_type', 
        'GI_score',
        'GI_pval', 
        'query_SMF', 
        'array_SMF',
        'double_mutant_fitness', 
        'double_mutant_fitness_std']


    # Load GI pairs
    df = pd.read_csv(fp, sep = '\t')
    df.columns = columns

    df = df.assign(query_orf=df.query_strain_ID.map(id_to_name))
    df = df.assign(array_orf=df.array_strain_ID.map(id_to_name))
    return df

import re
def get_best_strain(strain_gi_counts):
    strains = []
    counts = []
    subtypes = []
    prog = re.compile(r'[\w\d]+(-\w){0,1}_([a-zA-Z]+)\d+')
    # query_ids = arr_df.query_strain_ID.values
    # subtype = []

    for strain, count in strain_gi_counts.items():
        strains.append(strain)
        counts.append(count)
        m = prog.match(strain)
        subtypes.append( m.group(2))
    
    counts = np.asarray(counts)
    subtypes = np.asarray(subtypes)
    strains = np.asarray(strains)
    
    # split on count
    max_count_ind = counts == np.max(counts)
    if sum(max_count_ind) == 1:
        return strains[np.argmax(counts)]

    raise NotImplementedError('No tiebreaking rule implemented for strains: %s' % str(strain_gi_counts))

def get_best_strains(strain_ids, strain_orfs):
    orf_counts = {}
    # Get dictionaries 
    for strain, orf in tqdm(zip(strain_ids, strain_orfs)):
        if orf in orf_counts:
            if strain in orf_counts[orf]:
                orf_counts[orf][strain] += 1
            else:
                orf_counts[orf][strain] = 1
        else:
            orf_counts[orf] = {strain: 1}

    strains = [get_best_strain(counts) for orf, counts in orf_counts.items()]
    return strains


@click.command()
@click.option('-i', '--input', 'fp')
@click.option('--output_fmt')
@click.option('--GIs_only', is_flag=True)
def cli(fp, output_fmt, gis_only):
    # Renamed columns with shorter names...
    value_cols = ['GI_score',
                    'GI_pval', 
                    'query_SMF', 
                    'array_SMF',
                    'double_mutant_fitness', 
                    'double_mutant_fitness_std']
    if gis_only:
        # option for development, parse only GI scores and not all the other values
        value_cols = ['GI_score']

    # Load GI pairs
    _df = load_gis(fp)

    print(_df.head())
    print('* Read {} interactions'.format(len(_df)))
    print('\t- Total query genes: ', len(set(_df.query_orf)))
    print('\t- Total array genes: ', len(set(_df.array_orf)))


    print('* For each gene, getting query strain with most interactions')
    array_types = _df.array_type.unique()
    for array_type in array_types:
        print('\t- Getting array type: ', array_type)
        df = _df[_df.array_type == array_type]

        print('\t- getting best query strains')
        strains = get_best_strains(df.query_strain_ID, df.query_orf)
        df = df[df.query_strain_ID.isin(strains)]

        print('\t- getting best array strains')
        strains = get_best_strains(df.array_strain_ID, df.array_orf)
        df = df[df.array_strain_ID.isin(strains)]

        assert(len(set(df.query_strain_ID))  == len(set(df.query_orf)))
        assert(len(set(df.array_strain_ID))  == len(set(df.array_orf)))

        print('\t- Total query genes: ', len(set(df.query_orf)))
        print('\t- Total array genes: ', len(set(df.array_orf)))

        df = df.pivot(index='query_orf',columns='array_orf')
        for value_name in value_cols:
            gi_mat = df[value_name]
            print('\t-Extracting {} values (array type: {})'.format(value_name, array_type))
            print('\t\t- GIs shape:', gi_mat.shape)

            gi_mat = gi_mat.apply(pd.to_numeric)

            print('\t\t- GIs shape:', gi_mat.shape)
            print('\t\t- Sparsity:', sparsity(gi_mat.values))
            print('\t\t- Interactions measured :', np.sum(~np.isnan(gi_mat.values)))
            print('\t\t- Interactions missing:', np.sum(np.isnan(gi_mat.values)))
            
            gi_data = GIData(values=gi_mat.values.astype(float),
                        rows=gi_mat.index.values.astype(str),
                        cols=gi_mat.columns.values.astype(str),
                        check_symmetric=False)
            fp = output_fmt.format(value_name, array_type)

            gi_data.save(fp)
            print('\t\t- Saved values to:', fp)


        


if __name__ == "__main__":
    cli()