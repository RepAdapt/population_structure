"""Start population structure pipeline.

Usage
-----
In the following command, "/scratch" is the output directory (or a parent of the OUTDIR directory)

  module load apptainer
  
  sif=$HOME/pop_struct/population-structure.sif
  
  jobfile=$(
    apptainer exec -B "$HOME,/scratch:/scratch" $sif \
      conda run -n pop_struct \
      python $HOME/pop_struct/PopStruct_00_start_pipeline.py --vcf VCF -o OUTDIR [-h] [--discrete SAMP_TO_POP] \
  )
  
  cd $(dirname $jobfile) && sbatch $jobfile
  
Help
----
In the following command, "/scratch" is the output directory (or a parent of the OUTDIR directory)

  module load apptainer

  sif=$HOME/pop_struct/population-structure.sif

  apptainer exec -B "$HOME,/scratch:/scratch" $sif \
    conda run -n pop_struct \
    python $HOME/pop_struct/PopStruct_00_start_pipeline.py -h

"""
from pythonimports import *  # github.com/brandonlind/pythonimports

import argparse

import PopStruct_01_filter_vcf


def parse_command():
    parser = argparse.ArgumentParser(
      description="""RepAdapt population structure pipeline.

    Usage:
    ------
    sif=$HOME/pop_struct/population-structure.sif
      
    jobfile=$(
        apptainer exec -B "$HOME,/scratch:/scratch" $sif \\
          conda run -n pop_struct \\
          python $HOME/pop_struct/PopStruct_00_start_pipeline.py \\
            --vcf VCF \\
            -o OUTDIR \\
            [--discrete SAMP_TO_POP] \\
    )
    
    cd $(dirname $jobfile) && sbatch $jobfile
    """,
      add_help=False,
      formatter_class=argparse.RawTextHelpFormatter
    )
    
    requiredNAMED = parser.add_argument_group('required arguments')

    # -------- REQUIRED ARGS --------
    requiredNAMED.add_argument("--vcf",
                        required=True,
                        dest="vcf",
                        help='''Path to an unfiltered VCF. Assumed to be gzip-formatted. (Note downstream scripts assume unzipped vcf)''')

    requiredNAMED.add_argument("-o", "--output_directory",
                        required=True,
                        dest="outdir",
                        help='''Path to the directory for pipeline output.''')

    # -------- OPTIONAL ARGS --------
    parser.add_argument('-h', '--help',
                        action='help',
                        default=argparse.SUPPRESS,
                        help='Show this help message and exit.\n')

    parser.add_argument("--discrete",
                        required=False,
                        dest='samp_to_pop',
                        default='not_discrete',
                        help='''If populations are discrete, use this flag and point to file that maps sample names (column 1)
to **numeric** population IDs (column 2). Column names are arbitrary but must be included. This will be used to run hierfstat.''')
    ## --

    # PARSE ARGUMENTS
    args = parser.parse_args()

    if args.vcf:
        if not op.exists(args.vcf):
            raise Exception(f'The vcf does not exist in the specified path: {args.vcf}')
        # if not op.exists(f'{vcf}.tbi'):
        #     raise Exception(f'The vcf index does not exist: {vcf}.tbi')
    if not op.exists(args.outdir):
        os.makedirs(args.outdir, exist_ok=True)

    if args.samp_to_pop != 'not_discrete':
        assert op.exists(args.samp_to_pop)

    pkldump(args, f'{args.outdir}/pipeline_args.pkl')

    return args


def main():
    # parse arguments
    args = parse_command()

    PopStruct_01_filter_vcf.main(
        args.vcf,
        args.outdir,
        args.samp_to_pop
    )

    pass

if __name__ == '__main__':
    main()
