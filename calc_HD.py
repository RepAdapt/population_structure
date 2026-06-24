"""Calculate data to use for HD plots.

Usage
-----
python calc_HD.py vcf outdir

Parameters
----------
vcf : str | Path
    path to input vcf. can be .gz or unzipped
outdir : str | Path
    path to directory for saving dataframe txt file.
    basename of `vcf` is used as basename of output txt file.

Notes
-----
Writes a tab-delimited dataframe to /outdir/vcf_basename.txt
"""
import os
import sys
import vcfpy
import pandas as pd
from scipy.stats import binom
from collections import Counter

def get_het_depths(vcf, tmpfile):
    """Get allele depth stats from heterozygote individuals for specific loci of interest.
    
    modified from https://datadryad.org/resource/doi:10.5061/dryad.cm08m
    
    This is basically the dryad's vcf_to_allele_depth() from vcf_to_depth.py - but fixes errors in that script.
    """
    reader = vcfpy.Reader.from_path(vcf)
    lines = []
    snp_counter = 0
    for snp in reader:
        snp_counter += 1
        if snp_counter % 1000 == 0:
            print('snps processed =', snp_counter)
        
        if not snp.is_snv() or len(snp.ALT) > 1 or len(snp.ALT[0].value) > 1:  # must be a snv, biallelic, and not multinucleotide
            continue
    
        allele_counts = Counter()
        depth_a_of_ind, depth_b_of_ind = [], []
        total_depth_a, total_depth_b = 0, 0
        gtd_samps = 0  # genotyped samp count - added for revision
        for samp in snp.calls:
            if samp.called is True:  # if sample has a non-missing genotype - excludes eg './.' 
                gtd_samps += 1
    
                if samp.is_het is True:
                    refdepth, altdepth = samp.data['AD']
                    depth_a_of_ind.append(refdepth)
                    depth_b_of_ind.append(altdepth) 
                    total_depth_a += refdepth
                    total_depth_b += altdepth
    
                genotype = samp.data['GT']
                assert len(genotype) == 3, f'Unexpected genotype length. Expected string length=3. {genotype = }'
                allele_counts[genotype[0]] += 1  # genotype[0] takes first allele of a gt with sep='/' as well as sep='|' (phased)
                allele_counts[genotype[-1]] += 1
            
        num_hets = len(depth_a_of_ind)
        num_samples = gtd_samps  # fixed for manuscript revision - now reflects true sample count
        sum_a = sum(depth_a_of_ind)
        sum_b = sum(depth_b_of_ind)
        AF = allele_counts['1'] / (
            allele_counts['0'] + allele_counts['1']
        )
        MAF = AF if AF <= 0.5 else (1 - AF)
        if sum_a + sum_b > 0:
            line = '\t'.join([str(x) for x in [snp.CHROM, snp.POS, snp.REF, snp.ALT[0].value, AF, MAF,
                                               sum_a, sum_b , 
                                               sum_a / (sum_a + sum_b),
                                               num_hets, num_samples,
                                               total_depth_a, total_depth_b, 
                                               total_depth_a / (total_depth_a + total_depth_b),
                                               total_depth_b / (total_depth_a + total_depth_b)]])
            lines.append(line)
    
    print('done processing. total snps processed =', snp_counter)
    
    with open(tmpfile, 'w') as o:
        header = '\t'.join(['CHROM', 'POS', 'REF', 'ALT', 'AF', 'MAF',
                    'HET_REF_AD', 'HET_ALT_AD', 
                    'HET_REF_AD-div-HET_DP',
                    'num_hets', 'num_samples',
                    'TOTAL_REF_AD','TOTAL_ALT_AD',
                    'TOTAL_REF_AD-div-TOTAL_DP','TOTAL_ALT_AD-div-TOTAL_DP'])
        o.write("%s\n" % header)
        o.write("%s" % '\n'.join(lines))
    
    print(f'temp output written to {tmpfile}')
    return tmpfile


def get_z_scores(tmpfile, outfile):
    """Use read ratio statistics to calculate z-scores."""
    df = pd.read_table(tmpfile)

    # modified from https://datadryad.org/resource/doi:10.5061/dryad.cm08m : HDplot_python.py
    #SUM READ COUNTS PER LOCUS
    df['HET-TOTAL_AD'] = df['HET_REF_AD'] + df['HET_ALT_AD']
    
    #CALCULATE HETEROZYGOSITY perc
    df['hetPerc'] = [df.loc[row,'num_hets']/df.loc[row,'num_samples'] for row in df.index]
    
    #CALCULATE EXPECTED STANDARD DEVIATION BASED ON BINOMIAL DISTRIBUTION
    df['std'] = binom(n = df['HET-TOTAL_AD'], p = .5).std()
    
    #CALCULATE Z-SCORE BASED ON STANDARD DEVIATION
    df['z'] = (df['HET_REF_AD'] - (df['HET-TOTAL_AD']/2))/ df['std']

    df.to_csv(outfile, sep='\t', index=False, header=True)
    print(f'output written to {outfile}')

    pass


if __name__ == '__main__':
    thisfile, vcf, outdir = sys.argv

    outfile = f'{outdir}/' + os.path.basename(vcf).split('.vcf')[0] + '.txt'
    tmpfile = f'{outdir}/' + os.path.basename(outfile + '.tmp')

    get_het_depths(vcf, tmpfile)

    get_z_scores(tmpfile, outfile)

    if os.path.exists(outfile):
        os.remove(tmpfile)
        print('removed temp file')

    pass
    
