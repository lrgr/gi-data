import click
import pandas as pd

@click.group()
def cli():
    pass

@cli.command()
@click.option('-i', '--input_xls')
@click.option('-o', '--output')
@click.option('-c', '--usecols', type=int)
def xls_to_csv(input_xls, output, usecols):
    df = pd.read_excel(input_xls, usecols=range(usecols))
    df.to_csv(output, sep='\t', index=False)

@cli.command()
@click.option('-i', '--input_xls')
def cat_xls(input_xls):
    df = pd.read_excel(input_xls)
    print(df)

@cli.command()
@click.option('-i', '--input_xls', 'tsv')
@click.option('--sc-sp')
@click.option('--sp-sc')
def s3_homs(tsv, sc_sp, sp_sc):
    _sc = 'sc_orf'
    _sp = 'sp_orf'
    df = pd.read_csv(tsv, sep='\t')
    df[_sp] = df[_sp].str.upper()
    df[_sc] = df[_sc].str.upper()
    sc_sp_df = df[[_sc, _sp]]
    sp_sc_df = df[[_sp, _sc]]

    sc_sp_df.to_csv(sc_sp, index=False, header=False)
    sp_sc_df.to_csv(sp_sc, index=False, header=False)

    assert(sc_sp_df.shape == (239, 2))
    assert(sp_sc_df.shape == (239, 2))


if __name__ == "__main__":
    cli()
