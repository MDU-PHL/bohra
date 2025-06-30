from bohra.launcher.Utils import CustomFormatter, _check_path, _check_size_file
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


def _compare_r1_r2(r1:str, r2:str) -> bool:
    """Check if the R1 and R2 files are paired-end reads."""
    
    if r1 == "" or r2 == "":
        LOGGER.critical("For snippy you must provide both R1 and R2 files.")
        return False
    if _check_path(r1) and _check_path(r2):
        r1_size = _check_size_file(r1)
        r2_size = _check_size_file(r2)
        
        return True
    else:
        LOGGER.error("One or both of the read files do not exist.")
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
        


def _setup_comparative_args(kwargs:dict, command:dict) -> dict:
    
    command['modules'].append(kwargs['comparative_tool'])
    if kwargs['comparative_tool'] != "mash":
        command['params'].append(f"--tree_input {kwargs['tree_input']}")
    
    if kwargs['phylo'] and kwargs['comparative_tool'] != "mash":
        command['modules'].append(kwargs['tree_builder'])
    
    if kwargs['comparative_tool'] == "snippy" and not _check_input_snippy(kwargs['input_file']):
        LOGGER.critical("You have selected snippy as your tool for snp detection. Currently this requires paired-end reads. Please check your input file and try again.")
        raise SystemExit
    # command = _accessory_params(kwargs=kwargs, command=command)

    return command