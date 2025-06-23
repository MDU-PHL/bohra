import logging
import pathlib
import subprocess
import click
import os
import json

class CustomFormatter(logging.Formatter):
    """Logging Formatter to add colors and count warning / errors"""

    grey = "\x1b[38;21m"
    yellow = "\x1b[33;21m"
    red = "\x1b[31;21m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    format = "[%(levelname)s:%(asctime)s] %(message)s"

    FORMATS = {
        logging.DEBUG: grey + format + reset,
        logging.INFO: grey + format + reset,
        logging.WARNING: yellow + format + reset,
        logging.ERROR: red + format + reset,
        logging.CRITICAL: bold_red + format + reset
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt,datefmt='%m/%d/%Y %I:%M:%S %p')
        return formatter.format(record)
    


cfg_file = f"{pathlib.Path(__file__).parent.parent.resolve() / 'bohra_defaults.json'}"
with open(cfg_file, 'r') as f:
    CFG = json.load(f)


def _check_path(path):
    """
    Check if path provided exists and is accessible
    :input - path to where reads/contigs are
    :output - boolean 
    """ 
    
    if pathlib.Path(path).exists():
        
        return True
    else:
        
        raise SystemExit


def _run_subprocess(self, cmd):

    
    p = subprocess.run(cmd, shell = True)
    return p


def _get_required_columns() -> dict:
    """
    Get the required columns for the input file
    :return: dict of required columns
    """
    return {
        "required":["Isolate"],
        "must_have": ["r1","r2","assembly"],
        "optional": ["species"]
    }


def _get_columns_list() -> list:
    """
    Get the columns list for the input file
    :return: list of columns
    """
    return [
        "Isolate",
        "r1",
        "r2",
        "assembly"
    ]


def _get_pipelines(pipeline:str) -> list:
    """
    Get the list of pipelines available
    :return: list of pipelines
    """
    pipelines = {
        
        
    }

def _resource_opt() -> list:

    resource_options = [
        {   
            "name":"cpus",
            "help":"Number of max CPU cores to run, will define how many rules are run at a time, if 0 then the avail cpus will be determined at time of launch",
            "default":0,
        },
        {
            "name":"workdir",
            "help":"The directory where Bohra will be run, default is current directory",
            "default":pathlib.Path.cwd().absolute(),
            "type":click.Path(exists=True)
        },
        {
            "name":"conda_path",
            "help":"The path to where your pre-installed conda bohra-envs are stored. This can be provided in your profiles settings as well - it assumes you have pre-configured all of your conda environments for each process run by bohra, this is an advanced setting. Please take care if you are changing it.",
            "default":pathlib.Path(os.getenv('CONDA_PREFIX', ''))
        },
        {
            "name":"keep",
            "help":"If you are rerunning bohra over an exisiting directory set --keep to 'Y' to archive report files - otherwise previous report files will be removed.",
            "default":'N',
            "type":click.Choice(['Y', 'N'])
        },
        {
            "name":"force",
            "help":"Add if you would like to force a complete restart of the pipeline. All previous logs will be lost.",
            "is_flag":True,
            "default":False
        },
        {
            "name":"no_conda",
            "help":"Set if you DO NOT WANT to use separate conda environments for each nextflow process.",
            "is_flag":True,
            "default":False
        },
    
    ]

    return resource_options


def _get_common_options() -> list:
    """
    Get the common options available for each pipeline
    :return: list of options
    """
#     @click.option('--input_file', '-i',
#               help='Path to reads file, which is a tab-delimited with 3 columns <isolatename>  <path_to_read1> <path_to_read2>.',
#               default='')
# @click.option('--sylph_db', '-k',
#               default="",
#               show_default=True,
#               help="Path to DB for use with sylph")
# @click.option('--kraken_db', '-k',
#               default=os.getenv("KRAKEN2_DEFAULT_DB", ''),
#               metavar='KRAKEN2_DEFAULT_DB',
#               show_default=True,
#               help="Path to DB for use with kraken2")
# @click.option('--speciation/--no-speciation',
#               is_flag=True, 
#               help='Speciation will be performed by deafult - use --no-speciation if you do not need species detected.')

    common_options = [
        {
            "name":"input_file",
            "short_name":"-i",
            "help":"Path to tab-delimited input file. Required columns are 'Isolate', 'r1', 'r2', and 'assembly'. Optional column is 'species'. You can also supply additional columns of your choice for ttree annotation.",
            "default":'',
            "show_default":True
        },
        {
            "name":"sylph_db",
            "help":"Path to DB for use with sylph",
            "default":"",
            "show_default":True
        },
        {
            "name":"kraken_db",
            "help":"Path to DB for use with kraken2",
            "default":os.getenv("KRAKEN2_DEFAULT_DB", ''),
            "metavar":'KRAKEN2_DEFAULT_DB',
            "show_default":True
        },
        {
            "name":"speciation",
            "help":"Speciation will be performed by deafult. Use --no-speciation if you do not need species detected.",
            "is_flag":True,
            "default":True
        }
    ]

    return common_options


def _get_cmd_options() -> dict:
    """
    Get the click options available for each pipeline
    :return: dict of options
    """
    common_options = _get_common_options()
    # print(common_options)
    resource_options = _resource_opt()
    # print(resource_options)

    opt = common_options.copy()

    # opt.extend(resource_options)
    # print(opt)
    all_options = {
        "basic": [],
        "assemble": [
            {
                "name":"assembler",
                "short_name":"-a",
                "help":"Assembler to use (shovill uses spades > 3.14 < 4 with --isolate mode).",
                "default":"shovill",
                "type":click.Choice(['shovill', 'skesa', 'spades'])
            },
            {
                "name":"spades_args",
                "help":"Use to add arguments to spades (when running with --assembler spades) for example: '--cov-cutoff auto'",
                "default":""
            }
        ],
        "amr_typing": [
            {
                "name":"blast_db",
                "help":"Path to the mlst blast_db, defaults to what is installed in the environment.",
                "default":f"{os.getenv('BLAST_DB', '')}"
            },
            {
                "name":"data_dir",
                "help":"Path to the mlst datadir, defaults to what is installed in the environment.",
                "default":f"{os.getenv('PUBMLST_DB','')}"
            },
            {
                "name":"mlst_exclude",
                "short_name":"-me",
                "help":"mlst schemes to exclude - multiple possible ie -me scheme1 -me scheme2 -me scheme3",
                "default":[],
                "multiple":True
            },
            {
                "name":"mobsuite_db",
                "help":"Path to the mobsuite_db, defaults to what is installed in the bohra-mob_suite environment.",
                "default":f"{os.getenv('MOBSUITE_DB','')}"
            }
        ],
        "tb":[
            {
                "name":"reference_genome",
                "short_name":"-ref",
                "help":"The default reference genome for TB, if you want to use a different one please provide the path to the .fasta file.",
                "default":f"{pathlib.Path(__file__).parent.parent.resolve() / 'references' / 'tb' / 'NC_000962.3.fasta'}"
            },
            {
                "name":"mask",
                "short_name":"-m",
                "help":"Default mask file for TB, if you want to use a different one please provide the path to the .bed file.",
                "default":f"{pathlib.Path(__file__).parent.parent.resolve() / 'references' / 'tb' / 'mask.bed'}"
            },
            {
                "name":"cluster",
                "help":"Set if you want to do hierarchical clustering.",
                "is_flag":True,
                "default":True
            },
            {
                "name":"cluster_method",
                "short_name":"-cm",
                "help":"The clustering method to use, default is 'single-linkage'",
                "type":click.Choice(['single-linkage', 'average', 'complete', 'centroid', 'median', 'ward', 'weighted']),
                "default":"single-linkage"
            },
            {
                "name":"cluster_threshold",
                "short_name":"-ct",
                "help":"Comma separated list of thresholds to use for clustering",
                "type":str,
                "default":"5,12"
            },
        ],
        "full": [
            {   "name":"comparative_tool",
                "type":click.Choice(CFG["comparative_tools"]),
                "help": "Tool to use for comparative genomics.",
                "default": "snippy",
                "show_default": True
            },
            {
                "name":"reference_genome",
                "short_name":"-ref",
                "help":"Path to reference genome (.gbk or .fa). This is only required if using snippy for comparative analysis.",
                "default":""
            },
            {
                "name":"mask",
                "short_name":"-m",
                "help":"Path to mask file for snippy, if using snippy for comparative analysis.",
                "default":""
            },
            {
                "name":"gubbins",
                "help":"Set to use gubbins for recombination correction - only when using snippy.",
                "is_flag":True
            },
            {
                "name":"snippy_args",
                "help":"Additional arguments for snippy, if using snippy for comparative analysis.",
                "default":""
            },

            {
                "name":"minmap",
                "short_name":"-mp",
                "help":"Snippy - minimum read mapping quality to consider.",
                "default":"60"
            },
            {
                "name":"basequal",
                "short_name":"-bq",
                "help":"Snippy - Minimum base quality to consider.",
                "default":"13"
            },
            {
                "name":"minqual",
                "short_name":"-mq",
                "help":"Snippy - minimum QUALITY in VCF column 6",
                "default":"100"
            },
            {
                "name":"minfrac",
                "short_name":"-mf",
                "help":"Snippy - minimum proportion for variant evidence",
                "default":"0"
            },
            {
                "name":"ska_minfreq",
                "help":"Ska - minimum frequency for variant calling.",
                "default":"0.9"
            },
            {
                "name": "ska_alnargs",
                "help":"Ska - additional arguments for alignment.",
                "default":""
            },
            {
                "name":"ska2_kszise",
                "help":"Ska - kmer size for ska2, default is 31.",
                "default":"31"
            },
            {
                "name":"fuzzy_core_prop",
                "help":"Proportion of core genome to use for fuzzy core genome analysis.",
                "default":1.0,
                "type":float
            },
            {
                "name":"cluster",
                "help":"Set if you want to do hierarchical clustering.",
                "is_flag":True
            },
            {
                "name":"cluster_method",
                "short_name":"-cm",
                "help":"The clustering method to use, default is 'single-linkage'",
                "type":click.Choice(['single-linkage', 'average', 'complete', 'centroid', 'median', 'ward', 'weighted']),
                "default":"single-linkage"
            },
            {
                "name":"cluster_threshold",
                "short_name":"-ct",
                "help":"Comma separated list of thresholds to use for clustering, default is '10'",
                "type":str,
                "default":"10"
            },
            {
                "name":"no_phylo",
                "help":"Set if you do NOT want to generate a phylogenetic tree.",
                "is_flag":True
            },
            {
                "name":"tree_input",
                type:click.Choice(["distance", "alignment"]),
                "help":"Input type for tree building, either 'distance' or 'alignment'.",
                "default":"alignment"

            },
            {
                "name":"tree_builder",
                "help":"Tree builder to use for comparative analysis.",
                "default":"",
                "type":click.Choice(['veryfasttree', 'iqtree'])
            },
            
        ]
    }

    cmd_opt = {}
    for p in ["basic", "assemble", "amr_typing","full"]:  # Add other pipelines as needed
        opt.extend(all_options[p])
        popt = opt.copy()
        popt.extend(resource_options)
        cmd_opt[p] = popt
    
    comp = common_options.copy()
    comp.extend(all_options["full"])
    comp.extend(resource_options)
    cmd_opt["comparative"] = comp
    tb = common_options.copy()
    tb.extend(all_options["tb"])
    tb.extend(resource_options)
    cmd_opt["tb"] = tb

    return cmd_opt

# @click.option('--cpus',
#               help='Number of max CPU cores to run, will define how many rules are run at a time, if 0 then the avail cpus will be determined at time of launch', 
#               default=0)
# @click.option('--workdir', '-w',
#               default=pathlib.Path.cwd().absolute(), 
#               help='The directory where Bohra will be run, default is current directory', 
#               type=click.Path(exists=True))
# @click.option('--conda_path',       
#               default=pathlib.Path(os.getenv('CONDA_PREFIX', '')), 
#               help='The path to where your pre-installed conda bohra-envs are stored. This can be provided in your profiles settings as well - it assumes you have pre-configured all of your conda environments for each process run by bohra, this is an advanced setting. Please take care if you are changing it.')
# @click.option('--keep',
#               default='N', 
#               type=click.Choice(['Y', 'N']), 
#               help='If you are rerunning bohra over an exisiting directory set --keep to \'Y\' to archive report files - otherwise previous report files will be removed.')
# @click.option('--proceed',
#               is_flag=True, 
#               help='If you would like to proceed straigt to the pipeline.')
# @click.option('--force', '-f',
#               is_flag=True, 
#               help='Add if you would like to force a complete restart of the pipeline. All previous logs will be lost.')
# @click.option('--no-conda',
#               is_flag=True, 
#               help='Set if you DO NOT WANT to use separate conda environments for each nextflow process.')
# @click.option('--nfconfig','-nfcfg',
#               default = f"", 
#               help='An additional config file, required if running on a non-local machine, ie slurm, cloud. For help see documentation at https://github.com/MDU-PHL/bohra or https://www.nextflow.io/docs/latest/executor.html',) # don't need this
# @click.option('--profile',
#               default=f"", 
#               help='The resource profile to use. Defaults to local, if using an alternative config file, this value should represent the name of a profile provided')
