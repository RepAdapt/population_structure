"""Impute and LD-prune + clumping the SNPS.

TODO
----
- pass r2 thresh and window size from pipeline start pkl
"""

def main(vcf):
    basename = op.basename(vcf).split('.vcf')[0]

    text = f'''{sbatch_header}

hostname
date
    

echo PLINK

cd {impute_dir}

source $HOME/population_structure/conda_init.sh
conda activate lea_bigsnpr

plink --vcf {vcf} --make-bed --out {basename) --keep-allele-order --allow-extra-chr

date


echo IMPUTE

Rscript $HOME/population_structure/LEA_smnf_impute.R {vcf} {impute_dir}

date

'''

if __name__ == '__main__':
    pass