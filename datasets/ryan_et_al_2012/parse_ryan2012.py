import click
import pandas as pd
from GIData import GIData
from GIData.utils import sparsity

def map_names(map_dict, genes, verbose_indent='\t\t'):
    for i, gene in enumerate(genes):
        if gene in map_dict:
            new_name = map_dict[gene]

            print(verbose_indent, '- Renamed %s to %s' % (gene, new_name))
            genes[i] = new_name
    return genes

def has_prefix(prefixes, genes):
    prefixed_genes = []
    with_prefix = []

    for g in genes:
        for p in prefixes:
            if g.startswith(p) and not g == p:
                prefixed_genes.append(g)
                with_prefix.append(p)
    return list(zip(prefixed_genes, with_prefix))

@click.command()
@click.option('-i', '--input_emap')
@click.option('-o', '--output')
@click.option('-V', '--verbose', is_flag=True)

@click.option('--map_from', multiple=True, default=[])
@click.option('--map_to', multiple=True, default=[])


def parse_roguev(input_emap, output, verbose, map_from, map_to):
    df = pd.read_csv(input_emap, sep = '\t', header=0, index_col=0)
    df = df.astype(float)

    print('* Loading from: ', input_emap)
    rows = df.index.str.upper().values
    cols = df.columns.str.upper().values

    # check prefixes
    print('\t- genes in cols with prefix in rows')
    print(has_prefix(rows, cols))
    print('\t- genes in rows with predix in cols')
    print(has_prefix(cols, rows))

    assert len(map_from) == len(map_to)

    if len(map_from) > 0:
        map_dict = dict(zip(map_from, map_to))
        print('\t- Mapping (fixing) row genes')
        rows = map_names(map_dict, rows)
        print('\t- Mapping (fixing) col genes')
        cols = map_names(map_dict, cols)

    is_sym = len(rows) == len(cols) and all(rows == cols)
    print('\t- Data is symmetric?', is_sym)

    rowset = set(rows)
    colset = set(cols)
    l = rowset - colset
    r = colset - rowset
    n = rowset & colset
    print('\t- GIs shape:', df.shape)
    print('\t- rows : both : cols')
    print('\t- {} : {} : {}'.format(len(l), len(n), len(r)))
    if verbose:
        print('\t\t- row only genes:', sorted(l))
        print('\t\t- col only genes:', sorted(r))

    print('\t- Sparsity:', sparsity(df.values))
    print('\t- Rows have no duplicates?', len(rows) == len(rowset))
    print('\t- Cols have no duplicates?', len(cols) == len(colset))

    assert (len(rows) == len(rowset)), 'rows have duplicates'
    assert (len(cols) == len(colset)), 'cols have duplicates'

    gi_data = GIData(values=df.values.astype(float),
                     rows=df.index.values.astype(str),
                     cols=df.columns.values.astype(str),
                     check_symmetric=is_sym)
    
    print('* Saving to:', output)
    gi_data.save(output)

if __name__ == '__main__':
    parse_roguev()