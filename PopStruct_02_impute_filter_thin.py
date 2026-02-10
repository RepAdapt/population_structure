"""Impute, filter for MAF, and thin input VCF.

Input VCF assumed to be filtered for missing individuals and missing data.

Usage
-----
python PopStruct_02_impute_filter_thin.py vcf output_path is_discrete

Paramters
---------
vcf : Path
    path to uncompressed (.vcf) filtered vcf
output_path : Path
    path to output directory
is_discrete : 'not_discrete' | Path
    If populations are considered discrete, a path is specified to a file with the first column as the individual ID, and
    second column as the **numeric** population ID - this will run hierfstat in PopStruct_02_impute_filter_thin.py.
    Otherwise, hierfstat is not run - if not a path the default should be 'not_discrete'.

Notes
-----
the basename of the vcf file (eg basename.vcf) is used as a prefix for output files
"""
from pythonimports import *

def main(vcf, output_path, is_discrete, threads=32):
    if not vcf.endswith('.vcf'):
        raise Exception('Input VCF must be uncompressed (and end with .vcf)')

    sbatch_header = '\n'.join([line for line in read(f'%s/pop_struct/slurm_header_config.txt' % os.environ['HOME']) if line.startswith('#SBATCH')])
    
    basename = op.basename(vcf).split('.vcf')[0]

    job = f'{basename}_impute_filter_thin'

    hierfstat_text = ''
    if is_discrete != 'not_discrete':
        hierfstat_text += f'''echo HIERFSTAT
conda activate hierfstat

Rscript $HOME/code/GitHub/hierfstat.R {basename}_imputed_maf-filtered_thinned.txt {is_discrete}

date
'''

    text = f'''#!/bin/bash
#SBATCH --job-name={job}
#SBATCH --time=12:00:00
#SBATCH --mem=50000M
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={threads}
#SBATCH -o %x_%j.out
{sbatch_header}

hostname
date

echo PLINK

cd {output_path}

source $HOME/pop_struct/conda_init.sh
conda activate lea_bigsnpr

plink --vcf {vcf} --make-bed --out {basename} --keep-allele-order --allow-extra-chr

date


echo IMPUTE

Rscript $HOME/pop_struct/LEA_smnf_impute.R {vcf} {output_path} {threads}

date


echo LOSTRUCT

conda activate lostruct

Rscript $HOME/pop_struct/lostruct.R {basename}_imputed_maf-filtered.txt

date


{hierfstat_text}

'''
    shfile = f'{output_path}/shfiles/{job}.sh'
    with open(shfile, 'w') as o:
        o.write(text)
    sbatch(shfile)

    pass

if __name__ == '__main__':
    latest_commit(html=False)
    session_info.show(html=False)
    
    thisfile, vcf, output_path, is_discrete = sys.argv
    main(vcf, output_path, is_discrete)
    pass
