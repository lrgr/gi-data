import pandas as pd
import numpy as np
import click

import cloudpickle as cpkl
from scipy.stats import pearsonr
from sklearn.metrics import (r2_score, average_precision_score)

from GIData import GIData
from GIData.utils import sparsity

def cpkl_load(fp):
    with open(fp, 'rb') as f:
        return cpkl.load(f)

@click.command()
@click.option('--dcell_onto_scores')
@click.option('--dcell_onto_pvals')
@click.option('--costanzo10_scores')
@click.option('--costanzo10_pvals')
@click.option('--scores_output')
@click.option('--pvals_output')
@click.option('--mask_output')

def cli(dcell_onto_scores, dcell_onto_pvals, costanzo10_scores, costanzo10_pvals, 
        scores_output, pvals_output, mask_output):
    dcell_scores = cpkl_load(dcell_onto_scores)
    dcell_pvals = cpkl_load(dcell_onto_pvals)
    costanzo_scores = cpkl_load(costanzo10_scores)
    costanzo_pvals = cpkl_load(costanzo10_pvals)

    merged_scores, _ = merge_dcell_and_costanzo_data(dcell_scores, costanzo_scores)

    merged_pvals, dcell_only_mask = merge_dcell_and_costanzo_data(dcell_pvals, costanzo_pvals, 
                                                                  merged_rows= merged_scores.rows,
                                                                  merged_cols= merged_scores.cols)

    assert(np.all(merged_scores.cols == merged_pvals.cols))
    assert(np.all(merged_scores.rows == merged_pvals.rows))

    merged_scores.save(scores_output)
    merged_pvals.save(pvals_output)
    with open(mask_output, 'wb') as f:
        cpkl.dump(dict(values=dcell_only_mask, rows=merged_scores.rows, cols=merged_scores.cols), f)


def merge_dcell_and_costanzo_data(dcell_data, costanzo_data, merged_rows=None, merged_cols=None):
    
    dcell_rows, dcell_cols = dcell_data['rows'], dcell_data['cols']
    n_dcell_rows, n_dcell_cols = len(dcell_rows), len(dcell_cols)
    costanzo_rows, costanzo_cols = costanzo_data['rows'], costanzo_data['cols']

    if (merged_rows is None) and (merged_cols is None):
        # figure out rows and cols to add.
        new_rows_to_add = list(set(costanzo_rows) - set(dcell_rows))
        new_cols_to_add = list(set(costanzo_cols) - set(dcell_cols))

        merged_rows = np.append(dcell_rows, np.asarray(new_rows_to_add))
        merged_cols = np.append(dcell_cols, np.asarray(new_cols_to_add))

        assert(len(merged_rows) == len(costanzo_rows))
        assert(len(merged_cols) == len(costanzo_cols))
        assert(np.all(merged_rows[:n_dcell_rows] == dcell_rows))
        assert(np.all(merged_cols[:n_dcell_cols] == dcell_cols))
    elif (not merged_rows is None) and (not merged_cols is None):
        pass
    else:
        assert(False)

    # new reorder the costanzo data...
    costanzo_rows_g2i = dict((g,i) for i,g in enumerate(costanzo_rows))
    costanzo_cols_g2i = dict((g,i) for i,g in enumerate(costanzo_cols))

    row_reorder_idxs = [costanzo_rows_g2i[g] for g in merged_rows]
    col_reorder_idxs = [costanzo_cols_g2i[g] for g in merged_cols]

    merged_values = costanzo_data['values'].copy()
    merged_values = merged_values[row_reorder_idxs].copy()
    merged_values = merged_values[:, col_reorder_idxs].copy()

    # Insert dcell data into top left
    merged_values[: n_dcell_rows, : n_dcell_cols] = dcell_data['values']
    merged_data = GIData(values=merged_values.astype(float),
                         rows=merged_rows.astype(str),
                         cols=merged_cols.astype(str),
                         check_symmetric=False)

    mask = np.zeros_like(merged_values, dtype=bool)
    mask[: n_dcell_rows, : n_dcell_cols] = True

    return merged_data, mask

if __name__ == "__main__":
    cli()