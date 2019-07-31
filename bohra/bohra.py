"""Bohra 
.. moduleauthor:: Kristy Horan <kristyhoran15@gmail.com>

Bohra exinct tree kangaroo that lived on the nullarbor

Bohra is microbial genomics pipeline, designed predominantly for use in public health, but may also be useful in research settings. The pipeline takes as input a tab-delimited file with the isolate IDs followed by the path to READ1 and READ2, a reference for alignment and a unique identifier, where reads are illumina paired end reads (other platforms are not supported).

The pipline is based on nullarbor (https://github.com/tseemann/nullarbor) and is designed to be run in high performance computing environment.

Bohra is modular allowing the user to choose between calling SNPs and generating a phylogenetic tree, perform assemblies and detect AMR, perform typing etc; Or use the full pipeline to call SNPs, generate phylogenies, assemble, type and detect the pan-geneome. The output of Bohra is a html report that can be distributed, with downloadable tables and data.
"""


import logging
import argparse
import pathlib
import sys
import os
from bohra.SnpDetection import RunSnpDetection
from bohra.ReRunSnpDetection import ReRunSnpDetection



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


def main():
    # setup the parser
  
    parser = argparse.ArgumentParser(description='Bohra - a bacterial genomics pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # subparser for running the pipeline
    subparsers = parser.add_subparsers(help="Task to perform")

    parser_sub_run = subparsers.add_parser('run', help='Initial run of Bohra', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # options for running
    parser_sub_run.add_argument('--input_file','-i',help='Input file = tab-delimited with 3 columns <isolatename>  <path_to_read1> <path_to_read2>', default='')
    parser_sub_run.add_argument('--job_id','-j',help='Job ID, will be the name of the output directory', default='')
    parser_sub_run.add_argument('--reference','-r',help='Path to reference (.gbk or .fa)', default = '')
    parser_sub_run.add_argument('--mask','-m',default = False, help='Path to mask file if used (.bed)')
    parser_sub_run.add_argument('--pipeline','-p', default = 'sa', choices=['sa','s','a', 'all'], help=f"The pipeline to run. SNPS ('s') will call SNPs and generate phylogeny, ASSEMBLIES ('a') will generate assemblies and perform mlst and species identification using kraken2, SNPs and ASSEMBLIES ('sa' - default) will perform SNPs and ASSEMBLIES. ALL ('all') will perform SNPS, ASSEMBLIES and ROARY for pan-genome analysis")
    parser_sub_run.add_argument('--assembler','-a', default = 'shovill', choices=['shovill','skesa','spades'], help=f"Assembler to use.")
    parser_sub_run.add_argument('--cpus','-c',help='Number of CPU cores to run, will define how many rules are run at a time', default=36)
    parser_sub_run.add_argument('--minaln','-ma',help='Minimum percent alignment', default=0)
    parser_sub_run.add_argument('--prefillpath','-pf',help='Path to existing assemblies - in the form path_to_somewhere/isolatename/contigs.fa')
    parser_sub_run.add_argument('--mdu', default=True, help='If running on MDU data')
    parser_sub_run.add_argument('--workdir','-w', default = f"{pathlib.Path.cwd().absolute()}", help='The directory where Bohra will be run, default is current directory')
    parser_sub_run.add_argument('--resources','-s', default = f"{pathlib.Path(__file__).parent / 'templates'}", help='Directory where templates are stored')
    parser_sub_run.add_argument('--force','-f', action="store_true", help = "Add if you would like to force a complete restart of the pipeline. All previous logs will be lost.")
    parser_sub_run.add_argument('--dryrun','-n', action="store_true", help = "If you would like to see a dry run of commands to be executed.")
    parser_sub_run.add_argument('--gubbins','-g', action="store_true", help = "If you would like to run gubbins. NOT IN USE YET - PLEASE DO NOT USE")
    # parser for rerun
    
    parser_sub_rerun = subparsers.add_parser('rerun', help='Rerun of Bohra. Add or remove isolates from isolate list, change mask or reference.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # options for rerun
    parser_sub_rerun.add_argument('--reference','-r',help='Path to reference (.gbk or .fa)', default = '')
    parser_sub_rerun.add_argument('--mask','-m',default = '', help='Path to mask file if used (.bed)')
    parser_sub_rerun.add_argument('--cpus','-c',help='Number of CPU cores to run, will define how many rules are run at a time', default=36)
    parser_sub_rerun.add_argument('--workdir','-w', default = f"{pathlib.Path.cwd().absolute()}", help='Working directory, default is current directory')
    parser_sub_rerun.add_argument('--resources','-s', default = f"{pathlib.Path(__file__).parent / 'templates'}", help='Directory where templates are stored')
    parser_sub_rerun.add_argument('--dryrun','-n', action="store_true", help = "If you would like to see a dry run of commands to be executed.")
    parser_sub_rerun.add_argument('--gubbins','-g', action="store_true", help = "If you would like to run gubbins. NOT IN USE YET - PLEASE DO NOT USE")
    parser_sub_rerun.add_argument('--keep', '-k', action= 'store_true', help="Keep report from previous run")
    
    
    parser_sub_run.set_defaults(func=run_pipeline)
    
    parser_sub_rerun.set_defaults(func = rerun_pipeline)
    args = parser.parse_args()
    
    if vars(args) == {}:
        parser.print_help(sys.stderr)
    else:
        logging.basicConfig(filename='job.log',level=logging.INFO, format='[%(asctime)s] %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
        logger = logging.getLogger(__name__)
        args.func(args)
	
if __name__ == '__main__':
    main()

