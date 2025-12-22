from builtins import bool, dict, float, str
import logging
import pathlib
import subprocess
import click
import shutil
from shutil import SameFileError
import os
import json
import jinja2
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
    
# Logger
LOGGER =logging.getLogger(__name__) 
LOGGER.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(CustomFormatter())
fh = logging.FileHandler('bohra_test.log')
fh.setLevel(logging.DEBUG)
formatter = logging.Formatter('[%(levelname)s:%(asctime)s] %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p') 
fh.setFormatter(formatter)
LOGGER.addHandler(ch) 
LOGGER.addHandler(fh)

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




def _run_subprocess(cmd : str, capture_output: bool = True) -> subprocess.CompletedProcess:

    if capture_output:
        p = subprocess.run(cmd, shell = True, capture_output=True, encoding= 'utf-8')
        return p
    else:
        p = subprocess.run(cmd, shell = True)
        return p

def _get_required_columns() -> dict:
    """
    Get the required columns for the input file
    :return: dict of required columns
    """
    return {
        "required":["Isolate"],
        "must_have": ["r1","assembly"],
        "optional": ["r2", "species","is_control"]
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
        "assembly",
        "Species_expected",
        "is_control"
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
    pwd = f"{os.getenv('BOHRA_WD', pathlib.Path.cwd())}"
    resource_options = [
        {   
            "name":"cpus",
            "help":"Number of max CPU cores to run, will define how many rules are run at a time, if 0 then the avail cpus will be determined at time of launch",
            "default":1,
        },
        {
            "name":"workdir",
            "help":"The directory where Bohra will be run, default is current directory",
            "default":pwd,
            "type":click.Path(exists=True)
        },
        {
            "name":"use_conda",
            "help":"Use separate conda environments for each nextflow process.",
            "is_flag":True,
            "default":True
        },
        {
            "name":"dependency_prefix",
            "help":"The path to where your pre-installed conda bohra-envs are stored. This can be provided in your profiles settings as well - it assumes you have pre-configured all of your conda environments for each process run by bohra, this is an advanced setting. Please take care if you are changing it.",
            "default":f"{pathlib.Path( os.getenv('CONDA_PREFIX')) /  'bohra_conda_envs'}" if os.getenv('CONDA_PREFIX') else f"",
        },
        {
            "name":"profile_config",
            "help":"Path to the profile config file.",
            "default":f"{pathlib.Path( os.getenv('NF_PROFILE_CONFIG'))}" if os.getenv('NF_PROFILE_CONFIG') else f"",
        },
        
        {
            "name":"report_outdir",
            "short_name":"-ro",
            "help":f"The directory where Bohra will output results, default is report",
            "default":f"report",
            "help":"Please supply an output directory for bohra results.",
        },
        {
            "name":"replace-report",
            "help":"If you are rerunning bohra over an exisiting directory set --replace-report to override report files.",
            "is_flag":True,
            "default":False
        },
        {
            "name":"force",
            "help":"Add if you would like to force a complete restart of the pipeline. All previous logs and outputs will be lost.",
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

    prefix = f"{pathlib.Path(os.getenv('CONDA_PREFIX', ''))}"
    prefix = f"{pathlib.Path(prefix).name}" if prefix else 'bohra'

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
            "name":"kraken2_db",
            "help":"Path to DB for use with kraken2",
            "default":f"{pathlib.Path(os.getenv('KRAKEN2_DEFAULT_DB'))}" if os.getenv("KRAKEN2_DEFAULT_DB") else  "",
            "metavar":'BOHRA_KRAKEN2_DB',
            "show_default":True
        },
        {
            "name":"speciation",
            "help":"Speciation will be performed by deafult. Use none if you do not need species detected.",
            "type":click.Choice(['kraken2', 'none']),
            "default":"kraken2"
        },
        {
            "name":"no-auto-run",
            "default":False,
            "is_flag":True,
            "help":"Set --no-auto-run to prevent the pipeline from running automatically. You will need to copy and paste the command to run the pipeline yourself."
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
        },
        {
            "name":"trim/--no-trim",
            "help":"Set if you want to trim reads using fastp. Please note that this will duplicate reads. You may notice that you use alot more disk space if you are trimming.",
            "is_flag":True,
            "default":False
        },
        {
            "name":"no-downloadable-tables",
            "help":"Disable downloadable tables in the report html.",
            "is_flag":True,
            "default":False
        }
        ]
    

    return common_options


def _extract_tool_list(config_file:str, tool:str = "all")->dict:
    """
    Extract tool list from config file.
    """
    try:
        with open(config_file, 'r') as f:
            config = json.load(f)
            if tool != "all":
                return {tool: config[tool]}
            return config
    except Exception as e:
        LOGGER.critical(f"Error reading config file: {e}")
        raise SystemExit
    
def _get_dep_cmd_options() -> dict:
    """
    Get the click options available for each dependency
    :return: dict of options
    """
    envs_list = _extract_tool_list(f"{pathlib.Path(__file__).parent.parent}/config/dependencies.json")
    envs_list = list(envs_list.keys())
    deps_opts = {
        "install": [
            {
                "name": "tool",
                "type":click.Choice(envs_list + ["all"], case_sensitive=False),
                "default":"all",
                "help":"Install only a specific set of tools from a single environment. Should really only be used for development and/or testing purposes."
            }
        ],
        "update":[
            {
                "name": "tool",
                "type":click.Choice(envs_list + ["all"], case_sensitive=False),
                "default":"all",
                "help":"Update only a specific set of tools from a single environment. Should really only be used for development and/or testing purposes."
            }
        ],
        "check":[
            {
                "name": "tool",
                "type":click.Choice(envs_list + ["all"], case_sensitive=False),
                "default":"all",
                "help":"Update only a specific set of tools from a single environment. Should really only be used for development and/or testing purposes."
            }
        ]

    }

    return deps_opts

def _get_run_cmd_options() -> dict:
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
                "help":"Assembler to use (shovill_spades uses spades > 3.14 < 4 with --isolate mode).",
                "default":"shovill_spades",
                "type":click.Choice(['shovill_spades','shovill_skesa', 'skesa', 'spades'])
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
                "default":f"{os.getenv('BOHRA_BLAST_DB', '')}"
            },
            {
                "name":"data_dir",
                "help":"Path to the mlst datadir, defaults to what is installed in the environment.",
                "default":f"{os.getenv('BOHRA_PUBMLST_DB','')}"
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
                "default":f"{os.getenv('BOHRA_MOBSUITE_DB','')}"
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
                "name":"ignore_warnings",
                "help":"Set to ignore warnings during pipeline. Please note that this may lead to unexpected results.",
                "is_flag":True,
                "default":False
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
                "name":"cluster/--no-cluster",
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
                "type":click.Choice(["distance", "alignment"]),
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
                "name":"ignore_warnings",
                "help":"Set to ignore warnings during pipeline. Please note that this may lead to unexpected results.",
                "is_flag":True,
                "default":False
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
                "name":"cluster/--no-cluster",
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
    ],
        "mash":[],
    }

    for arg in comparative_args[tool]:
        if kwargs[arg] != "":
            command['params'].append(f"--{arg} {kwargs[arg]}")
    if kwargs['cluster']:
        command['params'].append(f"--cluster true")
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



def _get_template(template:str) -> jinja2.Template:
    """
    Function to get the template.
    Args:
        template (str): Path to the template file.
    Returns:
        jinja2.Template: Jinja2 template object.
    """
    if _check_path(template):
        with open(template, "r", encoding="utf-8") as file:
            return jinja2.Template(file.read())
    raise FileNotFoundError(f"Template file {template} not found.")


def _get_target(outpath:str, title:str) -> str:
    
    """
    Generate the target file path for an HTML file based on the given 
    output path and title.

    Args:
        outpath (str): The directory path where the file should be saved.
        title (str): The title of the file, which will be used to generate
        the filename.

    Returns:
        str: The full path to the target HTML file.

    Raises:
        FileNotFoundError: If the specified output path does not exist.
    """

    if pathlib.Path(outpath).exists():
        name = f"{title.replace(' ', '_').replace(':', '_').replace('/', '_').lower()}.nf.config"
        return pathlib.Path(outpath) / name
    raise FileNotFoundError(f"Output path {outpath} does not exist.")

