from bohra.launcher.Utils import CustomFormatter, _check_path, _check_size_file, _check_input_snippy
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


def _setup_mtb_args(kwargs:dict, command:dict, mtb:bool) -> dict:

    command['modules'].append('mtb')
    command['modules'].append(kwargs['comparative_tool'])
    command['params'].append(f"--tree_input alignment")
    
    if kwargs['phylo'] and kwargs['comparative_tool'] != "mash":
        command['modules'].append(kwargs['tree_builder'])
    if kwargs['comparative_tool'] == "snippy" and not _check_input_snippy(kwargs['input_file']):
        LOGGER.critical("You have selected snippy as your tool for snp detection. Currently this requires paired-end reads. Please check your input file and try again.")
        raise SystemExit



    return command

