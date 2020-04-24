import pandas as pd
import numpy as np
import click

import cloudpickle as cpkl
from scipy.stats import pearsonr
from sklearn.metrics import (r2_score, average_precision_score)

from GIData import GIData

def cpkl_load(fp):
    with open(fp, 'rb') as f:
        return cpkl.load(f)

from GIData.utils import sparsity

@click.command()
@click.option('--gis', 'data_fp', required=True)
@click.option('--pvals', 'pval_fp', required=False)
@click.option('--gis_output', required=True)
@click.option('--pvals_output', required=False)
def cli(data_fp, pval_fp, gis_output, pvals_output):
    print('* Loading data')

    print('\t- GIs from:', data_fp)

    gi_data = cpkl_load(data_fp)
    
    if pval_fp:
        print('\t- p-values from:', pval_fp)
        pval_data = cpkl_load(pval_fp)

    rows = gi_data['rows']
    cols = gi_data['cols']

    print('\t- # rows', len(rows))
    print('\t- # cols', len(cols))

    # all_genes = np.unique(np.concatenate((rows, cols)))
    # all_genes.sort()
    # N = len(all_genes)

    # print('\t- # total genes', N)

    row_g2i = dict((g,i) for i, g in enumerate(rows))
    col_g2i = dict((g,i) for i, g in enumerate(cols))

    asym_values = gi_data['values']
    values = asym_values.copy()
    print('\t- Sparsity: ', sparsity(asym_values))

    if pval_fp:
        asym_pvals = pval_data['values']
        pvals = asym_pvals.copy()

    for i, A in enumerate(rows):
        for j, B in enumerate(cols):
            # in row, col |  and col, row

            # indices in asym data
            A_r = row_g2i.get(A)
            B_c = col_g2i.get(B)
            A_c = col_g2i.get(A)
            B_r = row_g2i.get(B)

            if (A_r is not None) and \
               (B_c is not None) and \
               (A_c is not None) and \
               (B_r is not None):
                v_ij = asym_values[A_r, B_c]
                v_ji = asym_values[B_r, A_c]

                if v_ij * v_ji < 0:
                    # opposite signs...
                    values[A_r, B_c] = np.nan
                    values[B_r, A_c] = np.nan
                    if pval_fp:
                        pvals[A_r, B_c] = np.nan
                        pvals[B_r, A_c] = np.nan
                else:
                    if pval_fp:
                        p_ij = asym_pvals[A_r, B_c]
                        p_ji = asym_pvals[B_r, A_c]
                        p = min(p_ij, p_ji)
                        v = v_ij  if p_ij < p_ji else v_ji
                        pvals[A_r, B_c] = p
                        pvals[B_r, A_c] = p
                        values[A_r, B_c] = v
                        values[B_r, A_c] = v
                    else:
                        v = (v_ij + v_ji) / 2
                    
                        values[A_r, B_c] = v
                        values[B_r, A_c] = v
            # elif (A_r is not None) and \
            #      (B_c is not None):
            #     if pval_fp:
            #         p = asym_pvals[A_r, B_c]
            #     v = asym_values[A_r, B_c]
            # elif (A_c is not None) and \
            #      (B_r is not None):
            #     if pval_fp:
            #         p = asym_pvals[B_r, A_c]
            #     v = asym_values[B_r, A_c]
            # else:
            #     continue

            # values[i,j] = v
            # values[j,i] = v

            # if pval_fp:
            #     pvals[i,j] = p
            #     pvals[j,i] = p

    print('* Processed: summary stats')
    print('\t- Shape:', values.shape)
    print('\t- Sparsity: ', sparsity(values))
#    print('\t- total unique pairs:', np.sum(~np.isnan(values) / 2))

    processed_gis = GIData(values=values.astype(float),
                           rows=rows.astype(str),
                           cols=cols.astype(str),
                           check_symmetric=False)
    processed_gis.save(gis_output)

    if pval_fp:
        processed_pvals = GIData(values=pvals.astype(float),
                                 rows=rows.astype(str),
                                 cols=cols.astype(str),
                                 check_symmetric=False)
        processed_pvals.save(pvals_output)
    
if __name__ == "__main__":
    cli()