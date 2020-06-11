"""Bohra 
.. moduleauthor:: Kristy Horan <kristyhoran15@gmail.com>

Bohra exinct tree kangaroo that lived on the nullarbor

Bohra is microbial genomics pipeline, designed predominantly for use in public health, but may also be useful in research settings. The pipeline takes as input a tab-delimited file with the isolate IDs followed by the path to READ1 and READ2, a reference for alignment and a unique identifier, where reads are illumina paired end reads (other platforms are not supported).

The pipline is based on nullarbor (https://github.com/tseemann/nullarbor) and is designed to be run in high performance computing environment.

Bohra is modular allowing the user to choose between calling SNPs and generating a phylogenetic tree, perform assemblies and detect AMR, perform typing etc; Or use the full pipeline to call SNPs, generate phylogenies, assemble, type and detect the pan-geneome. The output of Bohra is a html report that can be distributed, with downloadable tables and data.
"""


import logging
import argparse
import configargparse
import pathlib
import sys
import os
import shutil
from bohra.SnpDetection import RunSnpDetection
from bohra.Utils import Nulla2bohra, UpdateBohra, CheckDeps
from bohra.version import version
# from SnpDetection import RunSnpDetection
# from Utils import Nulla2bohra, UpdateBohra, CheckDeps
# from version import version



#logging.basicConfig(filename='job.log',level=logging.INFO, format='[%(asctime)s] %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

def run_pipeline(args):
    '''
    Run the pipeline for the first time
    '''
    R = RunSnpDetection(args)
    return(R.run_pipeline())

def rerun_pipeline(args):
    '''
    Rerun the pipeline on a previous dataset, adding, removing isolates or changing reference or mask file
    '''
    R = ReRunSnpDetection(args)
    return(R.run_pipeline())

def nulla2bohra(args):
    '''
    ensure that bohra can be rerun over the top of an existing nullarbor directory
    '''
    N = Nulla2bohra(args)
    return(N.update())

def check_deps(args):
    """
    check that dependencies are all installed
    """
    C = CheckDeps(args)
    return(C.check())
    
def main():
    # setup the parser
  
    parser = configargparse.ArgumentParser(description=f'Bohra - a bacterial genomics pipeline - version {version}',formatter_class=configargparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + version)
    # subparser for running the pipeline
    subparsers = parser.add_subparsers(help="Task to perform")

    

    parser_sub_run = subparsers.add_parser('run', help='Initial run of Bohra', formatter_class=configargparse.ArgumentDefaultsHelpFormatter,default_config_files=[f"{pathlib.Path.cwd().absolute() / 'bohra.conf'}"])
    
    # options for running
    parser_sub_run.add_argument('--input_file','-i',help='Input file = tab-delimited with 3 columns <isolatename>  <path_to_read1> <path_to_read2>', default='')
    # parser_sub_run.add_argument('-S', '--use_singularity', action='store_true', help = 'Set if you would like to use singularity containers to run bohra.')
    # parser_sub_run.add_argument('--snippy_singularity', default='docker://mduphl/snippy:v4.4.5', help='The path to containers. If you want to use locally stored contianers please pull from dockerhub://mduphl/<toolname>.')
    # parser_sub_run.add_argument('--abritamr_singularity', default='docker://mduphl/abritamr:v0.2.2', help='The path to containers. If you want to use locally stored contianers please pull from dockerhub://mduphl/<toolname>.')
    parser_sub_run.add_argument('--job_id','-j',help='Job ID, will be the name of the output directory', default='')
    parser_sub_run.add_argument('--reference','-r',help='Path to reference (.gbk or .fa)', default = '')
    parser_sub_run.add_argument('--mask','-m',default = False, help='Path to mask file if used (.bed)')
    parser_sub_run.add_argument('--kraken_db', '-k', env_var="KRAKEN2_DEFAULT_DB", help="Path to DB for use with kraken2, if no DB present speciation will not be performed.")
    parser_sub_run.add_argument('--pipeline','-p', default = 'preview', choices=['preview','sa','all'], help=f"The pipeline to run. Preview (--preview - default) will calculate mash-distances and a mash-tree for quick inspection of your dataset. SNPs and ASSEMBLIES ('sa') will perform SNPs and ASSEMBLIES. ALL ('all') will perform SNPS, ASSEMBLIES and ROARY for pan-genome analysis")
    parser_sub_run.add_argument('--assembler','-a', default = 'shovill', choices=['shovill','skesa','spades'], help=f"Assembler to use.")
    parser_sub_run.add_argument('--cpus','-c',help='Number of CPU cores to run, will define how many rules are run at a time', default=16)
    parser_sub_run.add_argument('--minaln','-ma',help='Minimum percent alignment. Isolates which do not align to reference at this threshold will not be included in core phylogeny.', default=80)
    parser_sub_run.add_argument('--mincov','-mc',help='Minimum percent alignment. Isolates which do not have average read coverage above this threshold will not be included further analysis.', default=40)
    parser_sub_run.add_argument('--prefillpath','-pf',help='Path to existing assemblies - in the form path_to_somewhere/isolatename/contigs.fa')
    parser_sub_run.add_argument('-mdu', action = "store_true", help='If running on MDU data')
    parser_sub_run.add_argument('-workdir','-w', default = f"{pathlib.Path.cwd().absolute()}", help='The directory where Bohra will be run, default is current directory')
    parser_sub_run.add_argument('-resources','-s', default = f"{pathlib.Path(__file__).parent / 'templates'}", help='Directory where templates are stored')
    parser_sub_run.add_argument('-force','-f', action="store_true", help = "Add if you would like to force a complete restart of the pipeline. All previous logs will be lost.")
    parser_sub_run.add_argument('-dry-run','-n', action="store_true", help = "If you would like to see a dry run of commands to be executed.")
    parser_sub_run.add_argument('--cluster', action="store_true", help = "If you are running Bohra on a cluster.")
    parser_sub_run.add_argument('--json',help='Path to cluster.json - required if --cluster is set', default='')
    parser_sub_run.add_argument('--queue',help='Type of queue (sbatch or qsub currently supported) - required if --cluster is set.', default='')
    parser_sub_run.add_argument('--gubbins', '-g', action = 'store_true', help = 'Set to use gubbins for recombination correction.')
    parser_sub_run.add_argument('--keep', default = 'N', choices= ['Y', 'N'], help = 'If you are rerunning bohra over an exisiting directory set --keep to \'Y\' to archive report files - otherwise previous reprot files will be removed.')
    
    # parser for update
    parser_sub_nulla2bohra = subparsers.add_parser('nulla2bohra', help='Ensure that bohra can be rerun over an existing nullarbor folder. Can also be used to update older bohra directories. Must supply name of nullarbor directory, and your isolates.tab file', formatter_class=configargparse.ArgumentDefaultsHelpFormatter,default_config_files=[f"{pathlib.Path.cwd().absolute() / 'bohra.conf'}"])

    parser_sub_nulla2bohra.add_argument('-workdir','-w', default = f"{pathlib.Path.cwd().absolute()}", help='Working directory, default is current directory')
    parser_sub_nulla2bohra.add_argument('--job_id','-j',help='Job directory - the --name you used for nullarbor, will be the name of the output directory', default='')
    parser_sub_nulla2bohra.add_argument('--input_file','-i',help='Input file = tab-delimited with 3 columns <isolatename>  <path_to_read1> <path_to_read2>', default='')
    parser_sub_nulla2bohra.add_argument('--reference', '-r', help = 'The reference that was used in the previous run/ nullarbor job', default = '')
    
    # parser for checks
    parser_sub_check = subparsers.add_parser('check', help = "Check that all dependencies for bohra have been installed.")
    parser_sub_check.add_argument('--kraken_db', '-k', env_var="KRAKEN2_DEFAULT_DB", help="Path to DB for use with kraken2, if no DB present speciation will not be performed.")
    parser_sub_run.set_defaults(func=run_pipeline)
    parser_sub_nulla2bohra.set_defaults(func=nulla2bohra)
    parser_sub_check.set_defaults(func = check_deps)


        
    args = parser.parse_args()
    
    if len(sys.argv) <= 1:
        parser.print_help(sys.stderr)

    elif len(sys.argv) == 2 and sys.argv[1] in ['run', 'nulla2bohra']:
        parser.print_help(sys.stderr)
    
    else:
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.INFO)
        # create file handler which logs even debug messages
        fh = logging.FileHandler('bohra.log')
        fh.setLevel(logging.INFO)
        # create console handler with a higher log level
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        # create formatter and add it to the handlers
        formatter = logging.Formatter('[%(levelname)s:%(asctime)s] %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)
        # add the handlers to the logger
        logger.addHandler(ch)
        logger.addHandler(fh)
        logger.info(f"Starting bohra pipeline using {' '.join(sys.argv)}")
        args.func(args)
	
if __name__ == '__main__':
    main()

