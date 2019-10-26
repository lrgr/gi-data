import click

def load_homs(fp):
    '''Parse pombase homologs file'''
    homs = []
    import re
    with open(fp, 'r') as f:
        pat = re.compile(r'([^\(\)]+)(\([^\(\)]+\)|$)')
        for i, line in enumerate(f):
            elems= line.split('\t')
            if len(elems) != 2: continue
            sp_gene, sc_genes = elems
            sp_gene = sp_gene.upper()

            sc_genes = sc_genes.strip().split('|')
            sc_genes = [s.upper() for s in sc_genes]
            sc_genes = [pat.match(s).group(1) for s in sc_genes]
            for sc in sc_genes:
                if not sc == 'NONE':
                    homs.append((sp_gene, sc))
    return homs

@click.command()
@click.option('-i', '--input', 'fp')
@click.option('-o', '--output')
def cli(fp, output):
    homs = load_homs(fp)
    with open(output, 'w') as f:
        txt = '\n'.join(['%s\t%s' % pair for pair in homs])
        f.write(txt)

if __name__ == "__main__":
    cli()
