import click
import networkx as nx
import cloudpickle as cpkl
from itertools import product

def report_gis_in_homs(homs, left_gis, right_gis, left_name, right_name):
    combs = product(['rows', 'cols'], ['rows', 'cols'])
    msg = '{} ({}) x {} ({}) homologs:'
    frac_msg = '\tFraction of {} ({}): {:.3f}'
    for l, r in combs:
        hs = mwm_homs(homs, left_gis[l], right_gis[r])
        print(msg.format(left_name, l, right_name, r), len(hs))
        print(frac_msg.format(left_name, l, len(hs) / len(left_gis[l])))
        print(frac_msg.format(right_name, r, len(hs) / len(right_gis[r])))

def mwm_homs(homs, ls, rs):
    '''get max cardinamlity homs'''
    edgelist = []
    for l, r in homs:
        if l in ls and r in rs:
            edgelist.append((l,r))
    
    G = nx.from_edgelist(edgelist)
    _hs = nx.maximal_matching(G)
    hs = []
    for pair in _hs:
        if pair[0] in rs and pair[1] in ls:
            hs.append((pair[1], pair[0]))
        elif pair[1] in rs and pair[0] in ls:
            hs.append(pair)
        else:
            assert False, "We shouldn't have this case..."
    return hs

def load_gis(fp):
    with open(fp, 'rb') as f:
        return cpkl.load(f)

def load_homs(fp):
    homs = []
    with open(fp, 'r') as f:
        for line in f:
            homs.append(tuple(line.strip().split('\t')))
    l_homs, r_homs = zip(*homs)
    assert len(set(l_homs) & set(r_homs)) == 0
    return homs

def swap_cols(homs):
    l_homs, r_homs = zip(*homs)
    return list(zip(r_homs, l_homs))

def save_homs(homs, fp):
    print('\t- Saving homs to:', fp)
    with open(fp, 'w') as f:
        txt = '\n'.join(['%s\t%s' % pair for pair in homs])
        f.write(txt)

@click.command()
@click.option('-l', '--left_gis')
@click.option('-ln', '--left_name')
@click.option('-r', '--right_gis')
@click.option('-rn', '--right_name')
@click.option('-ho', '--homologs')
@click.option('-o', '--output_fmt')
def cli(left_gis, left_name, right_gis, right_name, homologs, output_fmt):
    l_gis = load_gis(left_gis)
    r_gis = load_gis(right_gis)
    combs = product(['rows', 'cols'], ['rows', 'cols'])
    
    homs = load_homs(homologs)

    print('* loaded {} homolog pairs (many-to-many)'.format(len(homs)))
    msg = '* {} ({}) x {} ({}) homologs:'
    frac_msg = '\t- Fraction of {} ({}): {:.3f}'
    for l_axis, r_axis in combs:
        hs = mwm_homs(homs, l_gis[l_axis], r_gis[r_axis])
        print(msg.format(left_name, l_axis, right_name, r_axis), len(hs))
        print(frac_msg.format(left_name, l_axis, len(hs) / len(l_gis[l_axis])))
        print(frac_msg.format(right_name, r_axis, len(hs) / len(r_gis[r_axis])))

        # Write to disk
        save_homs(hs, output_fmt.format(l=left_name, l_axis=l_axis, r=right_name, r_axis=r_axis))
        save_homs(swap_cols(hs), output_fmt.format(r=left_name, r_axis=l_axis, l=right_name, l_axis=r_axis))


if __name__ == "__main__":
    cli()
