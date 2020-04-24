import pandas as pd
import numpy as np
import click

import cloudpickle as cpkl
from scipy.stats import pearsonr
from sklearn.metrics import (r2_score, average_precision_score)

from GIData import GIData
from GIData.utils import sparsity

def sort_cols(df):
    # Swaps cols A and B if gene in B is < than gene in A
    swap = df['B'] < df['A']
    vals = df.loc[swap, ['B','A']].values
    df.loc[swap, ['A', 'B']] = vals
    return df

def report(df, key):
    print('On all GI scores')
    print('R2:', r2_score(df['true'], df[key]))
    print('pearsonr:', pearsonr(df['true'], df[key]))
    print('AUPR:', average_precision_score(df['true'] < -0.08, -1 * df[key]))
    print()
    print('On significant GI scores')
    sig_df = df[df['p'] < 0.05]
    print('R2:', r2_score(sig_df['true'], sig_df[key]))
    print('pearsonr:', pearsonr(sig_df['true'], sig_df[key]))
    print('AUPR:', average_precision_score(sig_df['true'] < -0.08, -1 * sig_df[key]))


@click.command()
@click.option('--dcell_fp')
@click.option('--onto_fp')
@click.option('--costanzo_2010_fp')
@click.option('--cpkl_output_fmt')
@click.option('--tsv_output')
def cli(onto_fp, dcell_fp, costanzo_2010_fp, cpkl_output_fmt, tsv_output):

    onto_raw_preds_fp = onto_fp
    df = pd.read_csv(onto_raw_preds_fp, sep='\t')
    df.columns = ['A', 'B', 'onto_pred', 'true', 'p']
    df = df[~df['true'].isna()].reset_index(drop=True)
    onto_df = df.sort_values(by=['A', 'B'])
    print(onto_df.head())
    print('Ontotype prediction data shape:', onto_df.shape)

    raw_dcell_preds = dcell_fp
    dcell_df = pd.read_csv(raw_dcell_preds, sep='\t', header=None)
    dcell_df.columns = ['A', 'B', 'dcell_pred']
    dcell_df = dcell_df.sort_values(by=['A', 'B'])
    print(dcell_df.head())
    print('DCell prediction data shape:', dcell_df.shape)
    onto_df = sort_cols(onto_df)
    dcell_df = sort_cols(dcell_df)

    keys = list(zip(onto_df['A'], onto_df['B']))
    keys = [a + ':' + b for (a,b) in keys]
    print('Ontotype preds have duplicate preds for gene pairs:',  not len(keys) == len(set(keys)))


    keys = list(zip(dcell_df['A'], dcell_df['B']))
    keys = [a + ':' + b for (a,b) in keys]
    print('DCell preds have duplicate preds for gene pairs:',  not len(keys) == len(set(keys)))

    print('Dropping duplicates in DCell preds')
    dcell_df = dcell_df.drop_duplicates(['A','B'])

    keys = list(zip(dcell_df['A'], dcell_df['B']))
    keys = [a + ':' + b for (a,b) in keys]
    print('DCell preds have duplicate preds for gene pairs:',  not len(keys) == len(set(keys)))

    print('Performing inner join to merge DCell and Ontotype predictions... (validating that join is 1-1)')
    merged_df = onto_df.merge(dcell_df, how='inner', left_on=['A', 'B'], right_on=['A', 'B'], validate='1:1')
    print(merged_df.head())
    print('** Ontotype: **')
    report(merged_df, 'onto_pred')


    print('** DCell: **')
    report(merged_df, 'dcell_pred')

    merged_df.to_csv(tsv_output, sep='\t', index=False)

    value_cols = ['onto_pred', 'dcell_pred', 'true', 'p']

    #TODO: we need to make this data match the shape of the data with our own processed Gs.
    # That is, we need to figure whether a pair is
    # A) query - array
    # b) array - query (so we swap)
    # c) or both array and query, then we add two values to the matrix.

    # check for duplicates
    # keys = list(zip(merged_df['A'], merged_df['B']))
    # keys = [frozenset((a,b)) for (a,b) in keys]
    # assert(len(keys) == len(set(keys)))
    # print('Merged scores are indexed by unique keys:', len(keys) == len(set(keys)))

    with open(costanzo_2010_fp, 'rb') as f:
        costanzo10_data = cpkl.load(f)
    
    costanzo10_rows = set(costanzo10_data['rows'])
    costanzo10_cols = set(costanzo10_data['cols'])
    print(len(costanzo10_rows), len(costanzo10_cols))

    needs_duplicates = df['A'].isin(costanzo10_rows) & df['A'].isin(costanzo10_cols) & \
                       df['B'].isin(costanzo10_rows) & df['B'].isin(costanzo10_cols)
    print(sum(needs_duplicates))

    # Incorrectly oriented (i.e. A is not in rows, or B is not in cols)
    needs_flipping = (~df['A'].isin(costanzo10_rows)) | (~df['B'].isin(costanzo10_cols))
    print(sum(needs_flipping))

    # Correcly oriented and does not need duplicating
    # keep =  df['A'].isin(costanzo10_rows) & (~df['A'].isin(costanzo10_cols)) & \
    #         df['B'].isin(costanzo10_cols) & (~df['B'].isin(costanzo10_rows))

    keep = ~(needs_flipping | needs_duplicates)
    print(sum(keep))

    dup_df = merged_df[needs_duplicates]
    dup_df.rename(columns={'A':'B', 'B':'A'}, inplace=True)

    flip_df = merged_df[needs_flipping]
    flip_df.rename(columns={'A':'B', 'B':'A'}, inplace=True)

    new_df = merged_df[needs_duplicates]
    new_df = new_df.append(dup_df)
    new_df = new_df.append(flip_df)
    new_df = new_df.append(merged_df[keep])
    print(new_df.head())
    print('row genes', len(set(new_df['A'])))
    print('col genes', len(set(new_df['B'])))

    print(len(new_df))

    new_df = new_df.pivot(index='A',columns='B')

    for value_name in value_cols:
        gi_mat = new_df[value_name]
        print('\t-Extracting {} values'.format(value_name))
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
        fp = cpkl_output_fmt.format(value_name)

        gi_data.save(fp)
        print('\t\t- Saved values to:', fp)

if __name__ == "__main__":
    cli()