"""Filter the VCF.

This will filter the VCF in the following order (not concurrently): 
    1) remove individuals with >10% missing data across SNPs, then
    2) remove SNPs with >10% missing data across individuals

Usage
-----
In the following command, "/scratch:/scratch" is an upper directory relative to OUTDIR
    module load apptainer
    
    sif=$HOME/pop_struct/population-structure.sif
    
    jobfile=$(
        apptainer exec -B "$HOME,/scratch:/scratch" $sif \
            conda run -n pop_struct \
            python $HOME/pop_struct/PopStruct_01_bcftools_filter.py vcf outdir is_discrete
    )
    
    cd $(dirname $jobfile) && sbatch $jobfile

Parameters
----------
vcf : Path
    Path to vcf file
outdir : Path
    Path to the output directory specified in PopStruct_00_start_pipeline.py
is_discrete : 'not_discrete' | Path
    If populations are considered discrete, a path is specified to a file with the first column as the individual ID, and
    second column as the **numeric** population ID - this will run hierfstat in PopStruct_02_impute_filter_thin.py. Column
    names are arbitrary but must be included in the file.
    If is_discrete is set to 'not_discrete' : hierfstat is not run/

TODO
----
- additional SBATCH flags - including email, partition, qos
- change to container instead of conda

"""
from pythonimports import *

def main(vcf, outdir, is_discrete):
    basename = op.basename(vcf).split('.vcf')[0]

    job = f'{basename}_filtered'

    shdir = makedir(f'{outdir}/shfiles')

    sbatch_header = '\n'.join([line for line in read(f'%s/pop_struct/slurm_header_config.txt' % os.environ['HOME']) if line.startswith('#SBATCH')])

    text = fr'''#!/bin/bash
#SBATCH --job-name={job}
#SBATCH --time=12:00:00
#SBATCH --mem=5000M
#SBATCH --ntasks=1
#SBATCH -o %x_%j.out
{sbatch_header}

hostname
date

module load apptainer

sif=$HOME/pop_struct/population-structure.sif

cd {outdir}

echo COMPUTE_IND_MISSINGNESS

apptainer exec -B "$HOME,/scratch:/scratch" $sif \
    conda run -n bcftools \
    bcftools stats -s - {vcf} \
    | awk '/^PSC/ {{nMissing=$14; total=$4+$5+$6+$14; miss=(total? nMissing/total : 0); if (miss<=0.10) print $3}}' \
    > keep.samples

date


echo FILTER_INDS

# keep only SNPs with a min and max number of alleles = 2
apptainer exec -B "$HOME,/scratch:/scratch" $sif \
    conda run -n bcftools \
    bcftools view -S keep.samples -m2 -M2 -v snps -Ou {vcf} \
    | bcftools +fill-tags -Ou -- -t F_MISSING \
    | bcftools view -i 'F_MISSING<=0.10' -Ou -o {job}.vcf

# tabix -p vcf {job}.vcf

date


echo SUBMIT_IMPUTATION

#conda activate pop_struct

#python $HOME/pop_struct/PopStruct_02_impute_filter_thin.py {job}.vcf {outdir} {is_discrete}

jobfile=$(
  apptainer exec -B "$HOME/pop_struct,/scratch:/scratch" $sif \
    conda run -n pop_struct \
    python $HOME/pop_struct/PopStruct_02_impute_filter_thin.py {job}.vcf {outdir} {is_discrete} | grep '\.sh$' | tail -n1
)

echo "jobfile=$jobfile"

date


echo SBATCHING

cd $(dirname $jobfile) && sbatch $jobfile

date

'''
    shfile = f'{shdir}/{job}.sh'
    with open(shfile, 'w') as o:
        o.write(text)
    # sbatch(shfile)
    print(shfile)

    pass



if __name__ == '__main__':
    latest_commit(html=False)
    session_info.show(html=False)

    thisfile, vcf, outdir, is_discrete = sys.argv

    pipeline_args = pklload(f'{outdir}/pipeline_args.pkl')

    main(vcf, outdir, is_discrete)
