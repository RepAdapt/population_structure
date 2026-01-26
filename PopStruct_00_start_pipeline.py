"""Start population structure pipeline.

TODO
----
- if is_filtered - bypass filtering
- if is_imputed no need to filter - bypass filtering and imputation

- config (or --not_slurm flag)
    - specify slurm is True or False - whether to sbatch the next step or just execute a command

"""
from pythonimports import *  # github.com/brandonlind/pythonimports

import argparse

import PopStruct_01_filter_vcf


def parse_command():
    parser = argparse.ArgumentParser(description="RepAdapt population structure pipeline.",
                                     add_help=False,
                                     formatter_class=argparse.RawTextHelpFormatter)
    
    requiredNAMED = parser.add_argument_group('required arguments')

    # REQUIRED ARGS
    requiredNAMED.add_argument("--vcf",
                        required=True,
                        dest="vcf",
                        help='''Path to a VCF. Assumed to be gzip-formatted. There should also be a bgzipped
format with the same basename but ending in .bgz (eg file.vcf.gz and file.vcf.bgz). The bgz file should be indexed
with tabix (tabix -p vcf file.vcf.bgz). If this VCF has been filtered for missingness across SNPs and individuals,
also use --filtered.''')

    requiredNAMED.add_argument("-o", "--output_directory",
                        required=True,
                        dest="outdir",
                        help='''Path to the directory for pipeline output.''')

    # -------- OPTIONAL ARGS --------
    parser.add_argument('-h', '--help',
                        action='help',
                        default=argparse.SUPPRESS,
                        help='Show this help message and exit.\n')

    ## -- determined if vcf file is filtered and imputed
    parser.add_argument('--filtered',
                        required=False,
                        action='store_false',
                        dest='is_filtered',
                        help='''Boolean: True if used (will not apply default filtering), False otherwise (will
apply default filtering). This specifies whether the VCF is filterd - bypasses filtering. Default False.''')
    parser.add_argument('--imputed',
                        required=False,
                        action='store_false',
                        dest='is_imputed',
                        help='''Boolean: True if used, False otherwise (impute with snmf() and impute() functions 
from the R package LEA v3.6.0. Whether the VCF has been imputed - bypasses filtering and imputation. Default False.''')
    ## --
    
    ## -- pruning / clumping flags --
    parser.add_argument('--r2_thresh',
                        required=False,
                        default=0.2,
                        dest='r2_thresh',
                        help='''R^2 threshold for pruning. Default R^2=0.2. See also --window_size
for important message.''')

    parser.add_argument('-w', '--window_size',
                        required=False,
                        default=500,
                        dest='window_size',
                        help='''The window size in BP. This is passed to the `size` argument of bigsnpr::snp_autoSVD.
Note - if --r2_thresh is not set to default, and this argument is not specified, 
bigsnpr `size` will default to 100 / r2_thresh. Default 500.''')
    ## --


    ## -- analysis args --

    parser.add_argument("--all",
                        required=False,
                        action='store_true',
                        dest='run_all',
                        help='''Boolean: True if used, False otherwise.
Whether to run all structure analyses; use this instead of using each analysis flag.''')
    ## --

    # in case running a structure analysis after others have been run, or applying different filtering
    parser.add_argument("--filtered_vcf",
                        required=False,
                        dest='filtered_vcf',
                        help='''Path to the pipeline-filtered VCF if one exists.''')

    # PARSE ARGUMENTS
    args = parser.parse_args()

    if args.vcf:
        if not op.exists(vcf):
            raise Exception(f'The vcf does not exist in the specified path: {args.vcf}')
        if not op.exists(f'{vcf}.tbi'):
            raise Exception(f'The vcf index does not exist: {vcf}.tbi')
    if not op.exists(args.outdir):
        os.makedirs(args.outdir, exist_ok=True)

    pkldump(args, f'{args.outdir}/pipeline_args.pkl')

    return args


def main():
    # parse arguments
    args = parse_command()

    # which structure analyses to run?
    if args.run_all is True:
        # set all analysis args to True
        pass

    if args.is_imputed is True:
        # bypass filtering and imputation
        PopStruct_02_impute.main(
            args.vcf,
            args.outdir
        )
        
    else:  # PS01 can bypass filtering
        PopStruct_01_filter_vcf.main(
            args.vcf,
            args.is_filtered,
            args.is_imputed,
            args.outdir
        )

    pass
    


if __name__ == '__main__':
    main()