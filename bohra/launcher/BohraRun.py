from bohra.launcher.Utils import CustomFormatter, _check_path, _run_subprocess, _get_required_columns, _get_columns_list
from bohra.launcher.BohraSetupFiles import _open_input_file, _check_data_format, _make_workdir
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

def _check_compute_resources(cpus:int,) -> bool:
    pass



def run_bohra(
        pipeline: str,
        kwargs: dict) -> bool:
    
    
    print(f"Running pipeline: {pipeline} with options: {kwargs}")

