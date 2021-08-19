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
    if args.check:
        return(R.check_dependencies())
    else:
        return(R.run_pipeline())

def check_deps(args):
    """
    check that dependencies are all installed
    """
    R = RunSnpDetection(args)
    return(R.check())
    
def main():
    # setup the parser
  
    parser = configargparse.ArgumentParser(description=f'Bohra - a bacterial genomics pipeline - version {version}',formatter_class=configargparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + version)
    # options for running
    parser.add_argument('--check',action="store_true", help = "Check that dependencies are installed correctly.")
    parser.add_argument('--input_file','-i',help='Path to reads file, which is a tab-delimited with 3 columns <isolatename>  <path_to_read1> <path_to_read2>. REQUIRED', default='')
    parser.add_argument('--contigs','-c',help='Path to contigs file, which is a tab-delimited with 3 columns <isolatename>  <path_to_contigs>. OPTIONAL if you already have assemblies.', default='')
    parser.add_argument('--job_id','-j',help='Job ID, will be the name of the output directory', default='')
    parser.add_argument('--reference','-r',help='Path to reference (.gbk or .fa)', default = '')
    parser.add_argument('--mask','-m',default = '', help='Path to mask file if used (.bed)')
    parser.add_argument("--abritamr_args",default="",help="Set if you would like to use point mutations, please provide a valid species.",       choices= ['Acinetobacter_baumannii', "Campylobacter", "Enterococcus_faecalis", "Enterococcus_faecium", "Escherichia", "Klebsiella", "Salmonella", "Staphylococcus_aureus", "Staphylococcus_pseudintermedius", "Streptococcus_agalactiae", "Streptococcus_pneumoniae", "Streptococcus_pyogenes", "Vibrio_cholerae"])
    parser.add_argument('--kraken_db', '-k', default="KRAKEN2_DEFAULT_DB", help="Path to DB for use with kraken2, if no DB present speciation will not be performed.")
    parser.add_argument('--pipeline','-p', default = 'preview', choices=['preview','default','all'], help=f"The pipeline to run. `preview` - generates a rapid tree using mash distances | `default` - runs snippy, phylogenetic tree (if > 3 sequences), assemblies, mlst and amr gene detection | `all` - same as default but includes roary pangenome analysis")
    parser.add_argument('--assembler','-a', default = 'spades', choices=['shovill','skesa','spades'], help=f"Assembler to use.")
    parser.add_argument('--cpus',help='Number of max CPU cores to run, will define how many rules are run at a time', default=16) # need to change
    parser.add_argument('--minaln','-ma',help='Minimum percent alignment. Isolates which do not align to reference at this threshold will not be included in core phylogeny.', default=0)
    parser.add_argument('--minqual','-mq',help='Minimum Avg quality of reads', default=0)
    parser.add_argument('--mincov','-mc',help='Minimum percent alignment. Isolates which do not have average read coverage above this threshold will not be included further analysis.', default=0)
    parser.add_argument('--workdir','-w', default = f"{pathlib.Path.cwd().absolute()}", help='The directory where Bohra will be run, default is current directory') # don't need this
    parser.add_argument('--force','-f', action="store_true", help = "Add if you would like to force a complete restart of the pipeline. All previous logs will be lost.")
    parser.add_argument('--no_phylo',action="store_true", help = "Set if you do NOT want to generate a phylogentic tree.")
    # parser.add_argument('--dry-run','-n', action="store_true", help = "If you would like to see a dry run of commands to be executed.")
    parser.add_argument('--executor',help='Type of queue (sbatch or qsub currently supported) - required if --cluster is set.', default='')
    # parser.add_argument('--gubbins', '-g', action = 'store_true', help = 'Set to use gubbins for recombination correction.')
    parser.add_argument('--keep', default = 'N', choices= ['Y', 'N'], help = 'If you are rerunning bohra over an exisiting directory set --keep to \'Y\' to archive report files - otherwise previous reprot files will be removed.')
    
    parser.set_defaults(func=run_pipeline)


        
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

