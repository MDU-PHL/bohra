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
from bohra.SnpDetection import RunSnpDetection
from bohra.ReRunSnpDetection import ReRunSnpDetection
from bohra.version import version


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
  
    parser = configargparse.ArgumentParser(description=f'Bohra - a bacterial genomics pipeline - version {version}',formatter_class=configargparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + version)
    # subparser for running the pipeline
    subparsers = parser.add_subparsers(help="Task to perform")

    

    parser_sub_run = subparsers.add_parser('run', help='Initial run of Bohra', formatter_class=configargparse.ArgumentDefaultsHelpFormatter,default_config_files=[f"{pathlib.Path.cwd().absolute() / 'bohra.conf'}"])
    
    # options for running
    parser_sub_run.add_argument('--input_file','-i',help='Input file = tab-delimited with 3 columns <isolatename>  <path_to_read1> <path_to_read2>', default='')
    parser_sub_run.add_argument('-S', '--use_singularity', action='store_true', help = 'Set if you would like to use singularity containers to run bohra.')
    parser_sub_run.add_argument('--singularity_path', default='shub://phgenomics-singularity', help='The path to singularity containers. If you want to use locally stored contianers please pull from shub://phgenomics-singularity (snippy.simg, prokka.simg, seqtk.simg, mash_kmc.simg, assemblers.simg, roary.simg). IMPORTANT bohra is designed to run with these containers... if you wish to use custom containers please contact developer or proceed at your own risk.')
    parser_sub_run.add_argument('--job_id','-j',help='Job ID, will be the name of the output directory', default='')
    parser_sub_run.add_argument('--reference','-r',help='Path to reference (.gbk or .fa)', default = '')
    parser_sub_run.add_argument('--mask','-m',default = False, help='Path to mask file if used (.bed)')
    parser_sub_run.add_argument('--kraken_db', '-k', env_var="KRAKEN2_DEFAULT_DB", help="Path to DB for use with kraken2, if no DB present speciation will not be performed.")
    parser_sub_run.add_argument('--pipeline','-p', default = 'sa', choices=['sa','s','a', 'all'], help=f"The pipeline to run. SNPS ('s') will call SNPs and generate phylogeny, ASSEMBLIES ('a') will generate assemblies and perform mlst and species identification using kraken2, SNPs and ASSEMBLIES ('sa' - default) will perform SNPs and ASSEMBLIES. ALL ('all') will perform SNPS, ASSEMBLIES and ROARY for pan-genome analysis")
    parser_sub_run.add_argument('--assembler','-a', default = 'shovill', choices=['shovill','skesa','spades'], help=f"Assembler to use.")
    parser_sub_run.add_argument('--cpus','-c',help='Number of CPU cores to run, will define how many rules are run at a time', default=36)
    parser_sub_run.add_argument('--minaln','-ma',help='Minimum percent alignment', default=0)
    parser_sub_run.add_argument('--prefillpath','-pf',help='Path to existing assemblies - in the form path_to_somewhere/isolatename/contigs.fa')
    parser_sub_run.add_argument('-mdu', action = "store_true", help='If running on MDU data')
    parser_sub_run.add_argument('-workdir','-w', default = f"{pathlib.Path.cwd().absolute()}", help='The directory where Bohra will be run, default is current directory')
    parser_sub_run.add_argument('-resources','-s', default = f"{pathlib.Path(__file__).parent / 'templates'}", help='Directory where templates are stored')
    parser_sub_run.add_argument('-force','-f', action="store_true", help = "Add if you would like to force a complete restart of the pipeline. All previous logs will be lost.")
    parser_sub_run.add_argument('-dry-run','-n', action="store_true", help = "If you would like to see a dry run of commands to be executed.")
    parser_sub_run.add_argument('--cluster', action="store_true", help = "If you are running Bohra on a cluster.")
    parser_sub_run.add_argument('--json',help='Path to cluster.json - required if --cluster is set', default='')
    parser_sub_run.add_argument('--queue',help='Type of queue (sbatch or qsub currently supported) - required if --cluster is set.', default='')
    
    # parser for rerun
    
    parser_sub_rerun = subparsers.add_parser('rerun', help='Rerun of Bohra. Add or remove isolates from isolate list, change mask or reference.', formatter_class=configargparse.ArgumentDefaultsHelpFormatter,default_config_files=[f"{pathlib.Path.cwd().absolute() / 'bohra.conf'}"])
    # options for rerun
    parser_sub_rerun.add_argument('-S', '--use_singularity', action='store_true', help = 'Set if you would like to use singularity containers to run bohra.')
    parser_sub_rerun.add_argument('--singularity_path', default='shub://phgenomics-singularity', help='The path to singularity containers. If you want to use locally stored contianers please pull from shub://phgenomics-singularity (snippy.simg, prokka.simg, seqtk.simg, mash_kmc.simg, assemblers.simg, roary.simg). IMPORTANT bohra is designed to run with these containers... if you wish to use custom containers please contact developer or proceed at your own risk.')
    parser_sub_rerun.add_argument('--reference','-r',help='Path to reference (.gbk or .fa)', default = '')
    parser_sub_rerun.add_argument('--mask','-m',default = '', help='Path to mask file if used (.bed)')
    parser_sub_rerun.add_argument('--cpus','-c',help='Number of CPU cores to run, will define how many rules are run at a time', default=36)
    parser_sub_rerun.add_argument('-workdir','-w', default = f"{pathlib.Path.cwd().absolute()}", help='Working directory, default is current directory')
    parser_sub_rerun.add_argument('--kraken_db', '-k', env_var="KRAKEN2_DEFAULT_DB", help="Path to DB for use with kraken2, if no DB present speciation will not be performed.")
    parser_sub_rerun.add_argument('-resources','-s', default = f"{pathlib.Path(__file__).parent / 'templates'}", help='Directory where templates are stored')
    parser_sub_rerun.add_argument('-dry-run','-n', action="store_true", help = "If you would like to see a dry run of commands to be executed.")
    parser_sub_rerun.add_argument('--gubbins','-g', action="store_true", help = "If you would like to run gubbins.")
    parser_sub_rerun.add_argument('-keep', action= 'store_true', help="Keep report from previous run")
    parser_sub_rerun.add_argument('-cluster', action="store_true", help = "If you are running Bohra on a cluster. Note if set you will need to provide a cluster.json file and a run_snakemake.sh, you can see examples on the documentation page.")
    parser_sub_rerun.add_argument('--json',help='Path to cluster.json - if not included will default to version provided in previous run', default='')
    parser_sub_rerun.add_argument('--queue',help='Type of queue (sbatch or qsub currently supported) - if not included will default to previous run', default='')
    
    parser_sub_run.set_defaults(func=run_pipeline)
    
    parser_sub_rerun.set_defaults(func = rerun_pipeline)
    args = parser.parse_args()
    
    if vars(args) == {}:
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

