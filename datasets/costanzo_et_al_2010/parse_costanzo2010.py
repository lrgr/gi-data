#!/usr/bin/env python

import numpy as np
import pandas as pd
import sys
import argparse

from os.path import join

from GIData import GIData
from GIData.utils import sparsity

#Parser
parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, required=True)
parser.add_argument('--output_fmt', type=str, required=True)

parser.add_argument('--with_cuttoffs', action='store_true')
parser.add_argument('--GIs_only', action='store_true')
parser.add_argument('--use_DAMP', action='store_true')
parser.add_argument('--use_ts_genes', action='store_true')
args=parser.parse_args()

if args.with_cuttoffs:
    '''
    The file contains the complete SGA genetic interaction dataset in a tab-delimited format with 7 columns:

    Query ORF
    Query gene name
    Array ORF
    Array gene name
    Genetic interaction score (ε)
    Standard deviation
    p-value
    '''
    columns = ['Query ORF',
                'Query gene name',
                'Array ORF',
                'Array gene name',
                'GI_score',
                'GI_std',
                'GI_pval']
    
    value_cols = ['GI_score', 'GI_std', 'GI_pval']
else:
    '''
    The file contains the complete SGA genetic interaction dataset in a tab-delimited format with 13 columns:

    Query ORF
    Query gene name
    Array ORF
    Array gene name
    Genetic interaction score (ε)
    Standard deviation
    p-value
    Query single mutant fitness (SMF)
    Query SMF standard deviation
    Array SMF
    Array SMF standard deviation
    Double mutant fitness
    Double mutant fitness standard deviation
    '''
    columns = ['Query ORF',
                'Query gene name',
                'Array ORF',
                'Array gene name',
                'GI_score',
                'GI_std',
                'GI_pval',
                'query_SMF',
                'query_SMF_std',
                'array_SMF', 
                'array_SMF_std',
                'double_mutant_fitness',
                'double_mutant_fitness_std']

    value_cols = ['GI_score',
                'GI_std',
                'GI_pval',
                'query_SMF',
                'query_SMF_std',
                'array_SMF', 
                'array_SMF_std',
                'double_mutant_fitness',
                'double_mutant_fitness_std']

if args.GIs_only:
    # option for development, parse only GI scores and not all the other values
    value_cols = ['GI_score']

# Load GI pairs
df = pd.read_csv(args.input,sep = '\t',header=None)
df.columns = columns

print('* Read {} interactions'.format(len(df)))

print()
print(df.head())
print()

# Convert to matrix
df = df.pivot(index='Query ORF',columns='Array ORF')

# Save each of the value types as seperate pickled matrices
for value_name in value_cols:
    gi_mat = df[value_name]
    print('* Extracting {}'.format(value_name))
    print('\t- GIs shape:', gi_mat.shape)

    # Remove DAMP genes
    if not args.use_DAMP:
        print('\t* Removing DAMP genes')
        damp_regex = '.+_damp'
        damp_cols = gi_mat.columns.str.match(damp_regex)
        damp_rows = gi_mat.index.str.match(damp_regex)
        print('\t\t- {} DAMP genes on cols'.format(np.sum(damp_cols)))
        print('\t\t- {} DAMP genes on rows'.format(np.sum(damp_rows)))
        gi_mat = gi_mat.iloc[~damp_rows, ~damp_cols]
        print('\t- GIs shape:', gi_mat.shape)
    
    # Remove TS alleles genes
    if not args.use_ts_genes:
        print('\t* Removing TS genes')
        ts_regex = '.+_tsq.*'
        ts_cols = gi_mat.columns.str.match(ts_regex)
        ts_rows = gi_mat.index.str.match(ts_regex)
        print('\t\t- {} TS genes on cols'.format(np.sum(ts_cols)))
        print('\t\t- {} TS genes on rows'.format(np.sum(ts_rows)))
        gi_mat = gi_mat.iloc[~ts_rows, ~ts_cols]
        print('\t- GIs shape:', gi_mat.shape)

    gi_mat = gi_mat.apply(pd.to_numeric)
    assert( len(gi_mat.index) == len(set(gi_mat.index)) )
    assert( len(gi_mat.columns) == len(set(gi_mat.columns)) )

    print('\t- GIs shape:', gi_mat.shape)
    print('\t- Sparsity:', sparsity(gi_mat.values))
    print('\t- Interactions measured :', np.sum(~np.isnan(gi_mat.values)))
    print('\t- Interactions missing:', np.sum(np.isnan(gi_mat.values)))
    
    gi_data = GIData(values=gi_mat.values.astype(float),
                  rows=gi_mat.index.values.astype(str),
                  cols=gi_mat.columns.values.astype(str),
                  check_symmetric=False)
    fp = args.output_fmt.format(value_name)
    gi_data.save(fp)
    print('\t- Saved values to:', fp)