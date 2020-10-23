''' generate manhattan plots '''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure


# functions to plot a simple manhattanplot which includes the bonferroni
# threshold and if applicable the permutation-based threshold


def manhattan(res_name, perm):
    ''' returns manhattan plot '''
    res = pd.read_csv(res_name).sort_values(['chr', 'pos'])
    if np.issubdtype(
            res['chr'].dtype,
            np.number) != True and np.issubdtype(
            res['chr'].dtype,
            np.number) != True:
        raise ValueError(
            '''The manhattan plot requires numeric
            information for the chromosomes and position of the markers
            ''')
    res.chr = res.chr.astype('category')
    res['my_cumsum'] = list(range(1, len(res) + 1))
    res['BPcum'] = 0
    my_s = 0
    bla = list()
    nbp = list()
    for i in res.chr.unique():
        nbp.append(np.max(res[res['chr'] == i]['pos']))
        bla.append(res[res['chr'] == i]['pos'] + s)
        my_s = my_s + nbp[i - 1]
    my_bla = [y for x in bla for y in x]
    res['BPcum'] = my_bla
    res['minuslog10pvalue'] = -np.log10(res.pval)
    res_group = res.groupby('chr')
    figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
    fig, my_axis = plt.subplots()
    del fig
    my_axis.margins(0.05)
    my_axis.hlines(
        np.log10(
            len(res)), xmin=0, xmax=np.max(
            res['BPcum']), linestyles="dashdot")
    if perm > 1:
        perm_res = pd.read_csv('perm_' + res_name)
        perm_idx = round(perm / 20, 0)
        perm_threshold = perm_res['min_p'][perm_idx]
        my_axis.hlines(-np.log10(perm_threshold), xmin=0,
                  xmax=np.max(res['BPcum']), linestyles="dashed")
    for name, group in res_group:
        my_axis.plot(
            group.BPcum,
            group.minuslog10pvalue,
            marker='o',
            linestyle='',
            ms=1,
            label=name)
    plt.xticks([])
    plt.legend()
    plt.savefig("manhattan.pdf")
