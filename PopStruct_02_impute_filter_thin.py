"""Impute, filter for MAF, and thin input VCF.

Input VCF assumed to be filtered for missing individuals and missing data.

Usage
-----
In the following command, "/scratch:/scratch" is an upper directory relative to outdir
    module load apptainer
    
    sif=$HOME/pop_struct/population-structure.sif
    
    jobfile=$(
        apptainer exec -B "$HOME/pop_struct,/scratch:/scratch" $sif \
            conda run -n pop_struct \
            python PopStruct_02_impute_filter_thin.py vcf outdir is_discrete
    )
    
    cd $(dirname $jobfile) && sbatch $jobfile
    
Paramters
---------
vcf : Path
    path to uncompressed (.vcf) filtered vcf
outdir : Path
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

def main(vcf, outdir, is_discrete, threads=32):
    outdir = op.abspath(outdir)  # in case it's a relative path

    if not vcf.endswith('.vcf'):
        raise Exception('Input VCF must be uncompressed (and end with .vcf)')

    sbatch_header = '\n'.join([line for line in read(f'%s/pop_struct/slurm_header_config.txt' % os.environ['HOME']) if line.startswith('#SBATCH')])
    
    basename = op.basename(vcf).split('.vcf')[0]

    job = f'{basename}_impute_filter_thin'

    hierfstat_text = ''
    if is_discrete != 'not_discrete':
        hierfstat_text += f'''echo HIERFSTAT
#conda activate hierfstat
#Rscript $HOME/pop_struct/hierfstat.R {basename}_imputed_maf-filtered.txt {is_discrete}
cd {outdir}
apptainer exec -B "$HOME,{outdir}:{outdir}" $sif \
    conda run -n hierfstat \
    Rscript $HOME/pop_struct/hierfstat.R {basename}_imputed_maf-filtered.txt {is_discrete}

date
'''

    text = fr'''#!/bin/bash
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

cd {outdir}

module load apptainer

sif=$HOME/pop_struct/population-structure.sif

apptainer exec -B "$HOME/pop_struct,{outdir}:{outdir}" $sif \
    conda run -n lea_bigsnpr \
    plink --vcf {vcf} --make-bed --out {basename} --keep-allele-order --allow-extra-chr --set-missing-var-ids @:# --double-id
date


echo IMPUTE

#Rscript $HOME/pop_struct/LEA_smnf_impute.R {vcf} {outdir} {threads}
apptainer exec -B "$HOME/pop_struct,{outdir}:{outdir}" $sif \
    env R_LIBS_USER="" R_LIBS="" \
    conda run --no-capture-output -n lea_bigsnpr \
    Rscript $HOME/pop_struct/LEA_smnf_impute.R {vcf} {outdir} {threads}

date


echo LOSTRUCT

#conda activate lostruct
#Rscript $HOME/pop_struct/lostruct.R {basename}_imputed_maf-filtered.txt

apptainer exec -B "$HOME/pop_struct,{outdir}:{outdir}" $sif \
    env R_LIBS_USER="" R_LIBS="" \
    conda run -n lea_bigsnpr \
    Rscript $HOME/pop_struct/lostruct.R {basename}_imputed_maf-filtered.txt

date


{hierfstat_text}

'''
    shfile = f'{outdir}/shfiles/{job}.sh'
    with open(shfile, 'w') as o:
        o.write(text)
    # sbatch(shfile)
    print(shfile)

    pass

if __name__ == '__main__':
    latest_commit(html=False)
    session_info.show(html=False)
    
    thisfile, vcf, outdir, is_discrete = sys.argv
    main(vcf, outdir, is_discrete)
    pass
