import click
import pandas as pd
from GIData import GIData
from GIData.utils import sparsity, multi_index_duplicated

import numpy as np

def get_orf(s):
    orf = s.split('(')[0]
    assert(orf[0] == 'S')
    assert(orf.find('.') >= 0)
    return s.split('(')[0]

@click.command()
@click.option('-i', '--input_emap')
@click.option('-o', '--output')
def parse_roguev(input_emap, output):
    '''
    In the Raw data Column and gene names look like:
    for deleted genes:
        - SPBC32F12.02(rec14SKI8)
    for DaMPed genes:
        - SPBC32F12.02(rec14SKI8) - DAMP
    The raw tsv file is symmetric, with the same rows and columns.
    One gene in the raw data is both DaMP and deleted. 
    We keep the deletion strain and drop the DAMPed strain.
    '''

    #Load initial GI matrix
    df = pd.read_csv(input_emap, sep = '\t', header=0, index_col=0)
    df = df.astype(float)

    #report initial matrix dimension
    print('Loaded initial %s x %s E-MAP' % df.shape)
    print('* Raw data:')
    print(df.head())
    print('\t- shape:', df.shape)

    assert(np.all(df.index.values == df.columns.values))

    # Array of DAMP/DEL labels
    mut_type_labels = np.repeat('DELETION', len(df))
    is_damped = df.index.str.endswith('DAMP')
    mut_type_labels[is_damped] = 'DAMP'

    genes = df.index.map(get_orf)
    genes = genes.str.upper()
    multi_index = pd.DataFrame(dict(Genes=genes, Mutation=mut_type_labels))
    multi_index = pd.MultiIndex.from_frame(multi_index,  names = ['Gene', 'Mutation'])
    df.columns = multi_index
    df.index = multi_index

    print('* Duplicated entries:', multi_index_duplicated(df.index, 'Gene'))
    del_df = df.loc[df.index.get_level_values('Mutation') == 'DELETION', 
                    df.columns.get_level_values('Mutation') == 'DELETION']
    
    print('* Keeping only deleted genes')

    gis = GIData.from_multiIndexDF(del_df)
    gis.save(output)
    
    print('* Processed GIs shape:', gis.shape)
    print('* Processed GIs sparsity:', sparsity(gis.values))

    pass

if __name__ == '__main__':
    parse_roguev()