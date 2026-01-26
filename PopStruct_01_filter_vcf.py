"""Filter the VCF.

This will filter the VCF in the following order (not concurrently): 
    1) remove individuals with >10% missing data across SNPs
    2) remove SNPs with >10% missing data across individuals

Usage: python PopStruct_01_bcftools_filter.py vcf is_filtered is_imputed outdir

Parameters
----------
vcf : Path
    Path to vcf file
is_filtered : bool
    Whether the vcf file has already been filtered
is_imputed : bool
    Whether the vcf file has already been imputed
outdir : Path
    Path to the output directory specified in PopStruct_00_start_pipeline.py

TODO
----
- FILL OUT LAST ELSE OF IFELSE block
- write command to file if not in slurm
- additional SBATCH flags - including email, partition, qos
- source conda_init

"""
from pythonimports import *

import PopStruct_02_impute


def main(vcf, is_filtered, is_imputed, outdir):
    basename = op.basename(vcf).split('.vcf')[0]

    filter_dir = makedir(f'{outdir}/00_filter_vcf')
    filter_outdir = makedir(f'{filter_dir}/filtered_vcf')
    shdir = makedir(f'{filter_dir}/shfiles')

    sbatch_header = f'''#!/bin/bash
#SBATCH --job-name={job}
#SBATCH --mem=5000M
#SBATCH --ntasks=1
#SBATCH -o %x_%j.out
'''

    if is_filtered is False:
        job = f'{basename}_filter'

        text = f'''{sbatch_header}
hostname
date

source {outdir}/conda_init.sh
conda activate bcftools

cd {filter_outdir}

echo COMPUTE_IND_MISSINGNESS

bcftools stats -s - {vcf} \
| awk '/^PSC/ {{nNonMissing=$5; nSites=$6; miss=(nSites-nNonMissing)/nSites; if(miss<=0.10) print $3}}' \
> keep.samples

date


echo FILTER_INDS

# keep only SNPs with a min and max number of alleles = 2
bcftools view -S keep.samples -m2 -M2 -v snps -Ou {vcf} \
| bcftools +fill-tags -Ou -- -t F_MISSING \
| bcftools view -i 'F_MISSING<=0.10' -Oz -o {job}.vcf.gz

tabix -p vcf {job}.vcf.gz

date


echo SUBMIT_IMPUTATION

python {pipeline_dir}/PopStruct_02_impute.py {filter_outdir}/{job}.vcf.gz

date

'''
        shfile = f'{shdir}/{job}.sh'
        with open(shfile, 'w') as o:
            o.write(text)
        sbatch(shfile)

    elif is_filtered is True and is_imputed is False:
        PopStruct_02_impute.main(vcf)

    else:
        #next step
        pass
    
    pass



if __name__ == '__main__':
    thisfile, vcf, is_filtered, is_imputed, outdir = sys.argv

    pipeline_dir = op.dirname(thisfile)

    pipeline_args = pklload(f'{outdir}/pipeline_args.pkl')

    main(vcf, is_filtered, is_imputed, outdir)