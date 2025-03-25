"""Bohra 
.. moduleauthor:: Kristy Horan <kristyhoran15@gmail.com>

Bohra exinct tree kangaroo that lived on the nullarbor

Bohra is microbial genomics pipeline, designed predominantly for use in public health, but may also be useful in research settings. The pipeline takes as input a tab-delimited file with the isolate IDs followed by the path to READ1 and READ2, a reference for alignment and a unique identifier, where reads are illumina paired end reads (other platforms are not supported).

The pipline is based on nullarbor (https://github.com/tseemann/nullarbor) and is designed to be run in high performance computing environment.

Bohra is modular allowing the user to choose between calling SNPs and generating a phylogenetic tree, perform assemblies and detect AMR, perform typing etc; Or use the full pipeline to call SNPs, generate phylogenies, assemble, type and detect the pan-geneome. The output of Bohra is a html report that can be distributed, with downloadable tables and data.
"""


# import logging
import argparse
import pathlib
import sys
import os
import shutil
from bohra.SnpDetection import RunSnpDetection, SetupInputFiles, InitBohra, TestBohra
from bohra.version import version


def run_pipeline(args):
    '''
    Run the pipeline for the first time
    '''
    
    R = RunSnpDetection(args)
    return(R.run_pipeline())

def find_reads(args):

    S = SetupInputFiles(args)
    S.find_reads()

def init_bohra():
    
    I = InitBohra()
    I.init()

def test_bohra(args):
    
    S = TestBohra(args)
    S.run_tests()
    # run_pipeline(inputs)

def main():
    # setup the parser
    # Check if CONDA_PREFIX is set
    conda_var = os.getenv('CONDA_PREFIX')
    if conda_var is None:
        # Set the default value
        conda_prefix = 'NO CONDA INSTALLED'
        # Set the default value
    else:
        conda_prefix = f"{pathlib.Path(os.getenv('CONDA_PREFIX')).parent}"
    
    parser = argparse.ArgumentParser(description=f'Bohra - a bacterial genomics pipeline - version {version}',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + version)
    # options for running
    subparsers = parser.add_subparsers(help="Task to perform")
    # run bohra pipeline
    parser_run = subparsers.add_parser('run', help='Run bohra', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser_run.add_argument('--check',action="store_true", help = "Check that dependencies are installed correctly.")
    parser_run.add_argument('--input_file','-i',help='Path to reads file, which is a tab-delimited with 3 columns <isolatename>  <path_to_read1> <path_to_read2>. REQUIRED', default='')
    parser_run.add_argument('--contigs','-c',help='Path to contigs file, which is a tab-delimited with 3 columns <isolatename>  <path_to_contigs>. OPTIONAL if you already have assemblies.', default='')
    parser_run.add_argument('--job_id','-j',help='Job ID is the name that will be displayed on your report', default='Bohra microbial genomics pipeline')
    parser_run.add_argument('--reference','-r',help='Path to reference (.gbk or .fa)', default = '')
    parser_run.add_argument('--mask','-m',default = '', help='Path to mask file if used (.bed)')
    parser_run.add_argument("--abritamr_args",default="",help="Set if you would like to use point mutations, please provide a valid species.", choices= ['Neisseria', 'Acinetobacter_baumannii', "Campylobacter", "Enterococcus_faecalis", "Enterococcus_faecium", "Escherichia", "Klebsiella", "Salmonella", "Staphylococcus_aureus", "Staphylococcus_pseudintermedius", "Streptococcus_agalactiae", "Streptococcus_pneumoniae", "Streptococcus_pyogenes", "Vibrio_cholerae"])
    parser_run.add_argument('--kraken_db', '-k', default="KRAKEN2_DEFAULT_DB", help="Path to DB for use with kraken2, if no DB present speciation will not be performed.")
    parser_run.add_argument('--pipeline','-p', default = 'preview', choices=['preview','default','full','snps','phylogeny','assemble','amr_typing'], help=f"The pipeline to run. `preview` - generates a rapid tree using mash distances | `default` - runs snippy, phylogenetic tree (if > 3 sequences), assemblies, mlst and amr gene detection | `all` - same as default but includes roary pangenome analysis")
    parser_run.add_argument('--assembler','-a', default = 'shovill', choices=['shovill','skesa','spades'], help=f"Assembler to use (shovill uses spades > 3.14 with --isolate mode).")
    parser_run.add_argument('--spades_args', default = "", help="Use to add arguments to spades (when running with --assembler spades) for example: '--cov-cutoff auto' ")
    parser_run.add_argument('--cpus',help='Number of max CPU cores to run, will define how many rules are run at a time, if 0 then the avail cpus will be determined at time of launch', default=0) # need to change
    parser_run.add_argument('--minmap','-mp',help='Snippy - minimum read mapping quality to consider.', default='60')
    parser_run.add_argument('--basequal','-bq',help='Snippy - Minimum base quality to consider.', default='13')
    parser_run.add_argument('--minqual','-mq',help='Snippy - minumum QUALITY in VCF column 6', default='100')
    parser_run.add_argument('--minfrac','-mf',help='Snippy - minumum proportion for variant evidence ', default='0')
    parser_run.add_argument('--mincov','-mc',help='Snippy - minimum site depth to for calling alleles.', default='10')
    parser_run.add_argument('--workdir','-w', default = f"{pathlib.Path.cwd().absolute()}", help='The directory where Bohra will be run, default is current directory') # don't need this
    parser_run.add_argument('--force','-f', action="store_true", help = "Add if you would like to force a complete restart of the pipeline. All previous logs will be lost.")
    parser_run.add_argument('--no_phylo',action="store_true", help = "Set if you do NOT want to generate a phylogentic tree.")
    parser_run.add_argument('--config', default = f"", help='An additional config file, required if running on a non-local machine, ie slurm, cloud. For help see documentation at https://github.com/MDU-PHL/bohra or https://www.nextflow.io/docs/latest/executor.html') # don't need this
    parser_run.add_argument('--profile', default = f"", help='The resource profile to use. Defaults to local, if using an alternative config file, this calue should represent the name of a profile provided') 
    parser_run.add_argument('--conda_path', default = conda_prefix, help='The path to where your pre-installed conda envs are stored, defaults to installing conda envs in your work directory. This can be provided in your profiles settings as well - it assumes you have pre-configured all of your conda environments for each process run by bohra, this is an advanced setting. Please take care.') 
    parser_run.add_argument('--blast_db', default = f"{os.getenv('BLAST_DB', '')}", help='Path to the mlst blast_db, defaults to what is installed in the environment.') 
    parser_run.add_argument('--data_dir', default = f"{os.getenv('PUBMLST_DB','')}", help='Path to the mlst datadir, defaults to what is installed in the environment.') 
    parser_run.add_argument('--mlst_exclude',default = f"", help='Space delimited list of mlst schemes to exclude.', nargs='+')
    parser_run.add_argument('--mobsuite_db', default = f"{os.getenv('MOBSUITE_DB','')}", help='Path to the mobsuite_db, defaults to what is installed in the environment.') 
    parser_run.add_argument('--proceed', action="store_true", help = "If you would like to proceed straigt to the pipeline.")
    parser_run.add_argument('--gubbins',  action = 'store_true', help = 'Set to use gubbins for recombination correction.')
    parser_run.add_argument('--keep', default = 'N', choices= ['Y', 'N'], help = 'If you are rerunning bohra over an exisiting directory set --keep to \'Y\' to archive report files - otherwise previous reprot files will be removed.')
    parser_run.add_argument('--no-conda',action="store_true", help="Set if you DO NOT WANT to use separate conda environments for each nextflow process.")

    parser_setup = subparsers.add_parser('generate_input', help='Generate input files for bohra', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_setup.add_argument(
        '--read_path',
        '-r',
        default='',
        help = 'Path to look for read files')
    parser_setup.add_argument(
        '--isolate_list',
        '-i',
        default= '',
        help = 'List of isolates (one isolate name per line) to include in input file. If not provided all sequences in found will be included.'
    )

    parser_init = subparsers.add_parser('init', help='Setup bohra dependencies', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_check = subparsers.add_parser('check', help='Check for bohra dependencies', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser_test = subparsers.add_parser('test', help='Test bohra', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_test.add_argument(
        '--reference',
        '-r',help='Reference to be used in test.', 
        default='Lm_Cluster1_J1-108.fa'
    )
    
    parser_run.set_defaults(func=run_pipeline)
    parser_setup.set_defaults(func=find_reads)
    parser_test.set_defaults(func=test_bohra)
    

        
    args = parser.parse_args()
    print(len(sys.argv))
    if len(sys.argv) < 2:
        parser.print_help(sys.stderr)
    elif sys.argv[1] in ['init', 'check']:
        init_bohra()
    else:
        args.func(args)
	
if __name__ == '__main__':
    main()

