import click
import pandas as pd
from GIData import GIData
from GIData.utils import sparsity

import numpy as np

@click.command()
@click.option('-i', '--input_emap')
@click.option('-o', '--output')

def parse_collins(input_emap, output):
    #Read raw data
    df = pd.read_csv(input_emap, sep ='\t',
                    header=None,
                    index_col=None)

    names = ['Gene', 'Mutation', 'Marker']

     # Get first three columns and first three rows for multi-index
    col_idxs = df.iloc[:3, 3:].T
    row_idxs = df.iloc[3:, :3]
    df = df.iloc[3:, 3:]
    col_multiIndex = pd.MultiIndex.from_frame(col_idxs, names = names)
    row_multiIndex = pd.MultiIndex.from_frame(row_idxs, names = names)
    df.index = row_multiIndex
    df.columns = col_multiIndex
    df = df.astype(np.float)

    # Some simple reporting
    print('* Processing data from: ', input_emap)
    print('* Raw data:')
    print(df.head())
    print('* Markers on columns:', df.columns.get_level_values('Marker').unique().values)
    print('* Markers on rows:', df.index.get_level_values('Marker').unique().values)
    print('* Mutations on columns:', df.columns.get_level_values('Mutation').unique().values)
    print('* Mutations on rows:', df.index.get_level_values('Mutation').unique().values)
    print()

    ### Get DELETIONS only dataframe:
    print('* Retrieving GI matrix with only deletions')

    del_df = df.loc[df.index.get_level_values('Mutation') == 'DELETION', 
                    df.columns.get_level_values('Mutation') == 'DELETION']
    
    gis = GIData.from_multiIndexDF(del_df)
    gis.save(output)
    
    print('* Processed GIs shape:', gis.shape)
    print('* Processed GIs sparsity:', sparsity(gis.values))

if __name__ == '__main__':
    parse_collins()