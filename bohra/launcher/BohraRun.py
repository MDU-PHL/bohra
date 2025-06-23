from bohra.launcher.Utils import CustomFormatter, _check_path, _run_subprocess, _get_required_columns, _get_columns_list
from bohra.launcher.BohraSetupFiles import _open_input_file, _check_data_format, _make_workdir
from bohra.launcher.BohraSetupResources import _get_profile, _set_cpu_limit_local
from bohra.launcher.BohraBasic import _setup_basic_args
from bohra.launcher.BohraAssembly import _setup_assembly_args
from bohra.launcher.BohraTyping import _setup_typing_args
import pandas as pd
import pathlib
import os
import logging

# Logger
LOGGER =logging.getLogger(__name__) 
LOGGER.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(CustomFormatter())
fh = logging.FileHandler('bohra_run.log')
fh.setLevel(logging.DEBUG)
formatter = logging.Formatter('[%(levelname)s:%(asctime)s] %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p') 
fh.setFormatter(formatter)
LOGGER.addHandler(ch) 
LOGGER.addHandler(fh)



def _setup_working_directory(input_file:str,
                             workdir : str) -> bool:

    input_data = _check_data_format(df = _open_input_file(input_file))

    _make_workdir(_input=input_data, workdir=workdir)
    
    
    return True

def _init_command_dict(profile:str, cpus:int) -> dict:

    return {"params":[
        f"--profile {profile}",
        f"-executor.cpus {cpus}",
        "-with-trace"
                ], "modules":[]}

def _check_assembly_req(kwargs: dict) -> bool:

    pass

def _funcs() -> dict:
    """Returns a dictionary of functions to be used in the pipeline."""
    
    return {
        "assemble": [_setup_assembly_args],
        "amr_typing":[ _setup_assembly_args, _setup_typing_args],
    }

def run_bohra(
        pipeline: str,
        kwargs: dict) -> bool:
    
    print(f"Running the {pipeline} pipeline with the following parameters:"
          f"\n{kwargs}")
    LOGGER.info(f"Checking on the setup for the {pipeline} pipeline.")
    profile = _get_profile(profile_config=kwargs.get('profile_config', ''),
                 profile=kwargs.get('profile', 'lcl'))
    max_cpus = int(_set_cpu_limit_local(cpus=kwargs.get('cpus', 0)))
    LOGGER.info(f"Using {int(max_cpus)} CPUs for the {pipeline} pipeline.")
    
    command = _init_command_dict(profile=profile, cpus=max_cpus)
    
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
        for _func in _funcs()[pipeline]:
            command = _func(kwargs=kwargs, command=command)
        # command = _funcs()[pipeline](kwargs=kwargs, command=command)

        print(f"Command to be run: {command}")
    else:
        LOGGER.error(f"Failed to create the working directory {kwargs['workdir']}.")
        raise SystemError
        