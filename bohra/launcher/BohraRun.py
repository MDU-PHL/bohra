from bohra.launcher.Utils import CustomFormatter, _check_path, _run_subprocess, _get_required_columns, _get_columns_list
from bohra.launcher.BohraSetupFiles import _open_input_file, _check_data_format, _make_workdir
from bohra.launcher.BohraSetupResources import _get_profile, _set_cpu_limit_local
from bohra.launcher.BohraBasic import _setup_basic_args
from bohra.launcher.BohraAssembly import _setup_assembly_args
from bohra.launcher.BohraTyping import _setup_typing_args
from bohra.launcher.BohraComparative import _setup_comparative_args
from bohra.launcher.BohraPangenome import _setup_pangenome_args
from bohra.launcher.CheckDeps import install_dependencies, check_databases

import pandas as pd
import pathlib
import os
import logging
import datetime
import subprocess

# Logger
LOGGER =logging.getLogger(__name__) 
LOGGER.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(CustomFormatter())
fh = logging.FileHandler('bohra.log')
fh.setLevel(logging.DEBUG)
formatter = logging.Formatter('[%(levelname)s:%(asctime)s] %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p') 
fh.setFormatter(formatter)
LOGGER.addHandler(ch) 
LOGGER.addHandler(fh)


def _check_deps_installed(prefix: str, check_only = False) -> bool:
    script_path = f"{pathlib.Path(__file__).parent}"
    
    deps = install_dependencies(prefix=prefix, check_only = check_only)
    if deps == 0:
        # LOGGER.info("All dependencies are installed.")
        return True
    else:
        LOGGER.error("One or more dependencies are not installed. Please run bohra check --install-deps to install them")
        raise SystemExit



def _check_keep(keep:str)->bool:

    if keep.lower() in ["yes", "y", "true", "t"]:
        now = datetime.datetime.today().strftime("%d_%m_%y_%H")
        try:
            cmd = f"cp -r {pathlib.Path.cwd() / 'report'} {pathlib.Path.cwd() / 'report_' + now + '_archive'}"
            proc = _run_subprocess(cmd=cmd)
            if proc.returncode == 0:
                LOGGER.info(f"Report directory archived successfully to {pathlib.Path.cwd() / 'report_' + now + '_archive'}.")
                return True
            else:
                LOGGER.error(f"Failed to archive the report directory. Command: {cmd}")
                return False
        except Exception as e:
            LOGGER.error(f"An error occurred while archiving the report directory: {e}")
            return False
    else:
        LOGGER.info("Report directory will be overridden by new results.")
        return True

def _setup_working_directory(input_file:str,
                             workdir : str) -> bool:

    input_data = _check_data_format(df = _open_input_file(input_file))

    _make_workdir(_input=input_data, workdir=workdir)
    
    
    return True

def _init_command_dict(profile:str, cpus:int, job_name:str, prefix:str, pipeline:str, profile_config:str, trim:bool) -> dict:

    command_dict = {"params":[
        f"-profile {profile}",
        f"--pipeline {pipeline}",
        f"-executor.cpus {cpus}",
        "-with-trace",
        f"--job_id {job_name}",
        f"--conda_prefix {prefix}",
                ], "modules":[]}

    if profile_config != "":
        command_dict["params"].append(f"-c {profile_config}")
    if trim:
        command_dict["modules"].append(f"trim")

    return command_dict

def _funcs() -> dict:
    """Returns a dictionary of functions to be used in the pipeline."""
    
    

    return {
        "basic":[],
        "assemble": [_setup_assembly_args],
        "amr_typing":[ _setup_typing_args],
        "full":[_setup_typing_args, _setup_comparative_args, _setup_pangenome_args ],
        "comparative":[ _setup_comparative_args ],
        "tb":[ _setup_comparative_args ],
    }

def _make_command(command: dict) -> str:
    """Constructs the command string from the command dictionary."""
    
    cmd = f"nextflow -Dnxf.pool.type=sync run {pathlib.Path(__file__).parent.parent.resolve() / 'bohra.nf'} " 
    
    cmd += " ".join(command["params"])
    if command["modules"]:
        cmd += " --modules " + ",".join(command["modules"])
    
    return cmd

def run_bohra(
        pipeline: str,
        kwargs: dict) -> bool:
    
    print(kwargs)
    LOGGER.info(f"Checking on the setup for the {pipeline} pipeline.")
    check_only = True if kwargs['install_deps'] == 'N' else False
    deps_ok = _check_deps_installed(prefix=kwargs.get('conda_prefix', 'bohra'), check_only = check_only)
    profile = _get_profile(profile_config=kwargs.get('profile_config', ''),
                 profile=kwargs.get('profile', 'lcl'))
    
    max_cpus = int(_set_cpu_limit_local(cpus=kwargs.get('cpus', 0)))
    LOGGER.info(f"Using {int(max_cpus)} CPUs for the {pipeline} pipeline.")

    command = _init_command_dict(profile=profile, profile_config = kwargs['profile_config'], cpus=max_cpus, job_name = kwargs.get('job_name', 'bohra'), prefix=kwargs.get('conda_prefix', 'bohra'), pipeline=pipeline, trim = kwargs["trim"])
    LOGGER.info(f"Checking if report directory needs to archived.")
    _check_keep(keep = kwargs["keep"])
    if _make_workdir(workdir=kwargs["workdir"],   
                  _input=kwargs["input_file"]):
        # update kwargs with checked input file
        kwargs["input_file"] = "input_checked.tsv"
    # add in the workdir and input file to the command to be run
        command["params"].append(f"--outdir {kwargs['workdir']}")
        LOGGER.info(f"Working directory {kwargs['workdir']} added to command successfully.")
        command["params"].append(f"--isolates {kwargs['input_file']}")
        LOGGER.info(f"Input file {kwargs['input_file']} added to command successfully.")
        command = _setup_basic_args(kwargs=kwargs, command=command)
        mtb = True if pipeline == "tb" else False
        LOGGER.info(f"Setting up arguments for the {pipeline} pipeline.")
        for _func in _funcs()[pipeline]:
            command = _func(kwargs=kwargs, command=command, mtb = mtb)
        # command = _funcs()[pipeline](kwargs=kwargs, command=command)
        
        cmd = _make_command(command=command)
        if kwargs["proceed"] == "N":
            LOGGER.info(f"Please paste the following command to run the pipeline:\n\033[1m{cmd}\033[0m")
        else:
            LOGGER.info(f"Running the command: {cmd}")
            proc = _run_subprocess(cmd=cmd)
            if proc.returncode == 0:
                LOGGER.info(f"The {pipeline} pipeline has completed successfully.")
                return True
            else:
                LOGGER.error(f"The {pipeline} pipeline failed with return code {proc.returncode}.")
                return False
    else:
        LOGGER.error(f"Failed to create the working directory {kwargs['workdir']}.")
        raise SystemError
        