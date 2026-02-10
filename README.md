# RepAdapt Population Structure pipeline
For any questions - contact Brandon Lind - lind.brandon.m@gmail.com

# Installation

1. Clone this repo on your file system: `git clone git@github.com:RepAdapt/population_structure.git`

    1. symlink the repo to your home directory: `ln -s /path/to/population_structure $HOME/pop_struct`

1. Install conda environments using a recent Anaconda version - eg Anaconda3-2025.12-1-Linux-x86_64.sh from the repo. 

    `bash Anaconda3-2025.12-1-Linux-x86_64.sh` - no need to initialize it in your ~/.bashrc.

    1. In the pipeline's `conda_init.sh` - edit and replace the path to your conda installation (likely `/your/home/anaconda3`).

    1. Source the pipeline's `conda_init.sh` to initialize Anaconda.

    1. Create needed environments.

        ```bash
        chmod +x conda_envs.sh
        ./conda_envs.sh
        ```

1. If your slurm system requires specialized #SBATCH flags (like `--partition` etc), add them to `slurm_header_config.txt`

1. If populations are discrete, create a file with two columns - the first is sample names, the second is a **numerical** population ID that will be used for hierfstat. In these cases use the --discrete flag to point to the file when starting the pipeline.

1. To start the pipeline (after sourcing `conda_init.sh`) run the following

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