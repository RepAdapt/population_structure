# RepAdapt Population Structure pipeline
For any questions - contact Brandon Lind - lind.brandon.m@gmail.com

# Outline

After starting the pipeline (see Usage), the following will be executed (submitted to slurm)
1. VCF filtering using bcftools (via a bash script created by PopStruct_01_filter_vcf.py)
    1. Samples are removed with > 10% missing data (threshold hard-coded)
    1. Then, loci are removed with > 10% missing data (threshold hard-coded)
1. Population structure is estimated (across several scripts, via a bash script created by PopStruct_02_impute_filter_thin.py)
    1. PLINK creates bim, bed, and fam files for locus and sample labelling within internal scripts that follow
    1. LEA_smnf_impute.R
        1. Dynamically choose K (explore K = 1-12) based on minimum cross-entropy via LEA::snmf using the filtered vcf in Step 1
            1. a K vs. CE plot is generated/saved for post hoc evaluation (suffix: ce_vs_k.pdf)
            1. an ancestry plot with the chosen K is generated/saved (suffix: ancestry_unimputed.pdf)
        1. Impute missing data using LEA::impute using best K from previous step
        1. Loci with MAF < 0.01 are removed
            1. A genotypes .txt file with samples for rownames and loci for column names is saved (suffix: imputed_maf-filtered.txt)
        1. bigsnpr::snp_autoSVD is used to LD prune loci (including long-range LD)
            1. A new K is dynamically chosen for this imputed/thinned SNP set
                1. a K vs. CE plot is generated saved for post hoc evaluation (suffix: imputed_maf-filtered_thinned_ce_vs_k.pdf)
            1. PCA coordinates are estimated and saved with rownames for samples and column names for PC (suffix: imputed_maf-filtered_thinned_pcs.txt)
            1. A genotypes .txt file with rownames for samples and column names is saved (suffix: imputed_maf-filtered_thinned.txt)
            1. Ancestry proportions are estimated/saved using the output from LD pruned loci
                1. Proportions saved (suffix: imputed_maf-filtered_thinned_ancestry_proportions.txt)
                1. A corresponding ancestry plot is generated (suffix: imputed_maf-filtered_thinned_ancestry.pdf)
    1. lostruct.R - using imputed + maf-filtered loci
        1. performs local PCA following directions on lostruct's [README](https://github.com/petrelharp/local_pca)
    1. hierfstat.R - using imputed + maf-filtered loci (if populations are designated as discrete)
        1. estimates/saves overall stats (suffix: hierfstat_overall_stats.txt)
        2. estimates/saves WC pairwise FST (suffix: hierfstat_pw_fst.txt)
        3. estimates/saves Ho He and Fis (suffix: hierfstat_mean_Ho_Hs_Fis.txt)

# Installation

1. Clone this repo on your file system: `git clone git@github.com:RepAdapt/population_structure.git`

    1. symlink the repo to your home directory: `ln -s /path/to/population_structure $HOME/pop_struct`

1. Install conda environments using a recent Anaconda version - eg Anaconda3-2025.12-1-Linux-x86_64.sh - you can get this by logging into Anaconda.com and clicking on Distribution. 

    `bash Anaconda3-2025.12-1-Linux-x86_64.sh` - no need to initialize it in your ~/.bashrc.

    1. In the pipeline's `conda_init.sh` - edit and replace the path to your conda installation (likely `/your/home/anaconda3`).

    1. Source the pipeline's `conda_init.sh` to initialize Anaconda.

    1. Create needed environments.

        ```bash
        chmod +x conda_envs.sh
        ./conda_envs.sh
        ```

1. If your slurm system requires specialized #SBATCH flags (like `--partition` etc), add them to `slurm_header_config.txt`

1. If populations are discrete, create a file with two columns - the first is sample names (same as in the unfiltered VCF), the second is a **numerical** population ID that will be used for hierfstat. In these cases use the --discrete flag to point to the file when starting the pipeline. Sample order does not matter - these are reorganized within the hierstat.R script to match genetic data. There are also several sanity checks (`stopifnot`) in the hierfstat.R script regarding sample set and order.

1. To start the pipeline (after sourcing `conda_init.sh` again) run the following

    ```bash
    conda activate pop_struct
    python PopStruct_00_start_pipeline.py --vcf VCF -o OUTDIR [-h] [--discrete SAMP_TO_POP]
    # see help menu for further info:
    python PopStruct_00_start_pipeline.py -h
    ```

    The pipeline will:

       1) filter outsamples and loci that have ≥ 10% missing data (see docstring at top of PopStruct_01_filter_vcf.py)

       2) Impute, filter for MAF, and thin input VCF (see docstring at the top of PopStruct_02_impute_filter_thin.py)

           - this calls on LEA_smnf_impute.R, hierfstat.R, lostruct.R

                - for more information and usage, see the docstrings of these files

                - to alter variables used for thinning, edit LEA_config.

     Note that even if not run on a slurm system, the script will still create the commands needed to run in the OUTDIR/shfiles directory.

1. At the end of each job, you can cat the slurm.out file for more information regarding outputs. See also Notes in docstring of hierfstat.R.

    `cat slurm_jobid.out | grep INFO`


# Usage

```
usage: PopStruct_00_start_pipeline.py --vcf VCF -o OUTDIR [-h] [--discrete SAMP_TO_POP]

RepAdapt population structure pipeline.

options:
  -h, --help            Show this help message and exit.
  --discrete SAMP_TO_POP
                        If populations are discrete, use this flag and point to file that maps sample names (column 1)
                        to **numeric** population IDs (column 2). Column names are arbitrary but must be included. This will be used to run hierfstat.

required arguments:
  --vcf VCF             Path to an unfiltered VCF. Assumed to be gzip-formatted. (Note downstream scripts assume unzipped vcf) 
  -o OUTDIR, --output_directory OUTDIR
                        Path to the directory for pipeline output.
```
