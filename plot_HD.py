"""Create figures using allele depth info from calc_HD.py.

Usage
-----
python plot_HD.py hd_file figdir [binsize]

Parameters
----------
hd_file : str | Path
    path to file created by calc_HD.py
figdir : str | Path
    path to save directory
[binsize] : int
    optional argument, used to adjust resolution of 2d scatter plots - default 100
    with few snps, lowering from 100 to around 40 or 50 may be appropriate
"""
import os
import sys
import myfigs as mf
import pandas as pd


def create_figs(hd_file, figdir, bins):
    """Create various figures, see titles.

    Parameters
    ----------
    hd_file : str | Path
        path to file created by calc_HD.py
    figdir : str | Path
        path to save directory
    """
    # print('reading hd_file')
    df = pd.read_table(hd_file)
    # print('done reading hd_file')
    basename = os.path.basename(hd_file).removesuffix('.txt')

    print('plotting read ratio deviation vs proportion of hets...')
    mf.scatter2d(
        x=df['hetPerc'],
        y=df['z'],
        ylab='Read Ratio Deviation\n(z-score values)',
        xlab='Proportion of Heterozygotes (H)',
        bins=bins,
        saveloc=f'{figdir}/{basename}_H-vs-Z.pdf',
        title='',
        title_kws={'fontsize' : 0}
    )

    print('plotting proportion of hets vs allelic ratio...')
    mf.scatter2d(
        x=df['hetPerc'],
        y=df['HET_REF_AD-div-HET_DP'],
        xlab='Proportion of heterozygotes (H)',
        ylab='Allelic Ratio (REF/total)',
        bins=bins,
        saveloc=f'{figdir}/{basename}_H-vs-AR.pdf',
        title='',
        title_kws={'fontsize' : 0}
    )

    print('plotting allelic ratio vs AF...')
    mf.scatter2d(
        x=df['AF'],
        y=df['HET_REF_AD-div-HET_DP'],
        ylab='Allelic Ratio (REF/total)',
        xlab='Frequency of ALT',
        bins=bins,
        saveloc=f'{figdir}/{basename}_AF-vs-AR.pdf',
        title='',
        title_kws={'fontsize' : 0}
    )

    print('plotting proportion of hets vs AF...')
    mf.scatter2d(
        x=df['AF'],
        y=df['hetPerc'],
        ylab='Proportion of heterozygotes (H)',
        xlab='Frequency of ALT',
        bins=bins,
        saveloc=f'{figdir}/{basename}_AF-vs-H.pdf',
        title='',
        title_kws={'fontsize' : 0}
    )

    pass


if __name__ == '__main__':
    thisfile, hd_file, figdir, *binsize = sys.argv

    if len(binsize) > 0:
        bins = binsize[0]
    else:
        bins = 100

    create_figs(hd_file, figdir, bins)

    pass
