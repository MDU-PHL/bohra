import logging
import pathlib
import subprocess
import click
import shutil
from shutil import SameFileError
import os
import json
import pandas as pd

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
    
    if pathlib.Path(path).exists() and os.access(path, os.R_OK):
        return True
    else:
        return False

def _check_size_file( path:str):
    '''
    check the size of a file
    '''
    s = path.stat().st_size
    return s




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
       "basic":"",
        "assemble":"assemble", 
        "amr_typing":"typing",
        "tb":"mtb"
        
    }

def _resource_opt() -> list:
    pwd = f"{pathlib.Path.cwd().absolute()}"
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
            "name":"use_conda",
            "help":"Use separate conda environments for each nextflow process.",
            "is_flag":True,
            "default":True
        },
        {
            "name":"conda_path",
            "help":"The path to where your pre-installed conda bohra-envs are stored. This can be provided in your profiles settings as well - it assumes you have pre-configured all of your conda environments for each process run by bohra, this is an advanced setting. Please take care if you are changing it.",
            "default":pathlib.Path(os.getenv('CONDA_PREFIX', ''))
        },
        {
            "name":"profile",
            "help":"The profile to use for running the pipeline. If not using defaults, you will need to have nextflow config set up - ist should be set as an environment variable - see docs for help. If not set, it will default to 'lcl'.",
            "default":os.getenv('TST_NF_PROFILE', 'lcl'),
        },
        {
            "name":"profile_config",
            "help":"Path to the profile config file.",
            "default":f"{pathlib.Path( os.getenv('TST_NFCFG', ''))}", # need to change this to reflect the ENV variable or the default supplied with bohra
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
        
    
    ]

    return resource_options


def _get_common_options() -> list:
    """
    Get the common options available for each pipeline
    :return: list of options
    """

    common_options = [
        {
            "name":"input_file",
            "short_name":"-i",
            "help":"Path to tab-delimited input file. Required columns are 'Isolate', 'r1', 'r2', and 'assembly'. Optional column is 'species'. You can also supply additional columns of your choice for ttree annotation.",
            "default":'',
            "show_default":True
        },
        {
            "name":"job_name",
            "short_name":"-j",
            "help":"Name of the job to be run, this will be used to name the output report html.",
            "default":"bohra",
        },
        {
            "name":"sylph_db",
            "help":"Path to DB for use with sylph",
            "default":os.getenv("SYLPH_DEFAULT_DB", ''),
            "show_default":True
        },
        {
            "name":"kraken2_db",
            "help":"Path to DB for use with kraken2",
            "default":os.getenv("KRAKEN2_DEFAULT_DB", ''),
            "metavar":'KRAKEN2_DEFAULT_DB',
            "show_default":True
        },
        {
            "name":"speciation",
            "help":"Speciation will be performed by deafult. Use none if you do not need species detected.",
            "type":click.Choice(['kraken2', 'sylph', 'none']),
            "default":"sylph"
        },
        {
            "name":"text_color",
            "help":"Color to use for the text in the report html. Default is 'white'.",
            "default":"#ffffff",
        },
        {
            "name":"background_color",
            "help":"Color to use for the background in the report html. Default is '#343a40'.",
            "default":"#343a40",
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
            {   "name":"comparative_tool",
                "type":click.Choice(CFG["comparative_tools"]),
                "help": "Tool to use for comparative genomics.",
                "default": "snippy",
                "show_default": True
            },
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
                "name":"snippy_args",
                "help":"Additional arguments for snippy, if using snippy for comparative analysis.",
                "default":""
            },

            {
                "name":"minmap",
                "short_name":"-mp",
                "help":"Snippy - minimum read mapping quality to consider.",
                "default":60
            },
            {
                "name":"basequal",
                "short_name":"-bq",
                "help":"Snippy - Minimum base quality to consider.",
                "default":13
            },
            {
                "name":"minqual",
                "short_name":"-mq",
                "help":"Snippy - minimum QUALITY in VCF column 6",
                "default":100
            },
            {
                "name":"minfrac",
                "short_name":"-mf",
                "help":"Snippy - minimum proportion for variant evidence",
                "default":0
            },
            {
                "name":"fuzzy_core_prop",
                "help":"Snippy - proportion of core genome to use for fuzzy core genome analysis.",
                "default":1.0,
                "type":float
            },
            {
                "name":"ska_minfreq",
                "help":"Ska - minimum frequency for variant calling.",
                "default":0.9
            },
            {
                "name": "ska_alnargs",
                "help":"Ska - additional arguments for alignment.",
                "default":""
            },
            {
                "name":"ska2_kszise",
                "help":"Ska - kmer size for ska2, default is 31.",
                "default":31
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
                "type":click.Choice(['single', 'average', 'complete', 'centroid', 'median', 'ward', 'weighted']),
                "default":"single",
    
            },
            {
                "name":"cluster_threshold",
                "short_name":"-ct",
                "help":"Comma separated list of thresholds to use for clustering",
                "type":str,
                "default":"5,12"
            },
            {
                "name":"phylo",
                "help":"Set if you do want to generate a phylogenetic tree.",
                "is_flag":True,
                "default":True
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
                "default":"veryfasttree",
                "type":click.Choice(['veryfasttree', 'iqtree']),
                
            },
            {
                "name":"annotations",
                "help":"Comma separated list of annotations to use for the tree, default is 'Tx:cluster_threshold'. MUST be present in the input file.",
                "default":"",
            }
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
                "is_flag":True,
                "default":False
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
                "default":60
            },
            {
                "name":"basequal",
                "short_name":"-bq",
                "help":"Snippy - Minimum base quality to consider.",
                "default":13
            },
            {
                "name":"minqual",
                "short_name":"-mq",
                "help":"Snippy - minimum QUALITY in VCF column 6",
                "default":100
            },
            {
                "name":"minfrac",
                "short_name":"-mf",
                "help":"Snippy - minimum proportion for variant evidence",
                "default":0
            },
            {
                "name":"fuzzy_core_prop",
                "help":"Snippy - proportion of core genome to use for fuzzy core genome analysis.",
                "default":1.0,
                "type":float
            },
            {
                "name":"ska_minfreq",
                "help":"Ska - minimum frequency for variant calling.",
                "default":0.9
            },
            {
                "name": "ska_alnargs",
                "help":"Ska - additional arguments for alignment.",
                "default":""
            },
            {
                "name":"ska2_kszise",
                "help":"Ska - kmer size for ska2, default is 31.",
                "default":31
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
                "help":"The clustering method to use, default is 'single'",
                "type":click.Choice(['single', 'average', 'complete', 'centroid', 'median', 'ward', 'weighted']),
                "default":"single"
            },
            {
                "name":"cluster_threshold",
                "short_name":"-ct",
                "help":"Comma separated list of thresholds to use for clustering, default is '10'",
                "type":str,
                "default":"10"
            },
            {
                "name":"phylo",
                "help":"Set if you do want to generate a phylogenetic tree.",
                "is_flag":True,
                "default":True
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
                "default":"veryfasttree",
                "type":click.Choice(['veryfasttree', 'iqtree']),
                
            },
            {
                "name":"annotations",
                "help":"Comma separated list of annotations to use for the tree, default is 'Tx:cluster_threshold'. MUST be present in the input file.",
                "default":"",
            },
            {
                "name":"pangenome_groups",
                "short_name":"-pg",
                "help":"Groups for analysis of pangenome data",
                "type":click.Choice(["clusters","mlst"]),
                "default":"clusters"
            }
            
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


def _compartive_args(tool:str,kwargs:dict, command:dict) -> dict:
    """
    Get the comparative arguments available for each pipeline
    :return: list of options
    """
    comparative_args = {
        "snippy":["snippy_args",
            "minmap",
            "basequal",
            "minqual",
            "minfrac",
            "fuzzy_core_prop"],
        "ska":["ska_minfreq",
               "ska_alnargs",
            "ska2_kszise"
    ]
    }

    for arg in comparative_args[tool]:
        if kwargs[arg] != "":
            command['params'].append(f"--{arg} {kwargs[arg]}")
    if kwargs['cluster']:
        for arg in ["cluster_method", "cluster_threshold"]:
            if kwargs[arg] != "":
                command['params'].append(f"--{arg} {kwargs[arg]}")
    return command

 

def _compare_r1_r2(r1:str, r2:str) -> bool:
    """Check if the R1 and R2 files are paired-end reads."""
    
    if r1 == "" or r2 == "":
        return False
    if _check_path(r1) and _check_path(r2):
        return True
    else:
        return False    

def _check_input_snippy(input_file:str) -> bool:
    """Check if the input file is a valid snippy output."""
    reads = False
    if _check_path(input_file):
        df = pd.read_csv(input_file, sep='\t')
        df = df.fillna('')
        for row in df.iterrows():
            chk = _compare_r1_r2(r1=row[1]['r1'], r2=row[1]['r2'])
            if  chk:
                reads = True
    
    return reads

def _check_reference(ref:str) -> bool:
    """Check if the reference genome is a valid file."""
    
    if ref != "" and _check_path(ref):
        try:
            dst = pathlib.Path(ref).name
            shutil.copy(ref, dst)
            ref = f"{pathlib.Path.cwd() / dst}"
        except SameFileError:
            ref = f"{pathlib.Path.cwd() / dst}"
        return ref
    else:
        return "no_ref"
    
def _check_mask(mask:str) -> bool:
    """Check if the mask is a valid file."""
    msk = "no_mask"
    if  mask != "" and _check_path(mask):
        try:
            dst = pathlib.Path(mask).name
            shutil.copy(mask, dst)
            msk = f"{pathlib.Path.cwd() / dst}"
        except SameFileError:
            msk = f"{pathlib.Path.cwd() / dst}"
        
        msk = f"{pathlib.Path.cwd() / dst}"
    
    return msk