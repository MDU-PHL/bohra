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




def _setup_assembly_args(kwargs:dict, command:dict) -> dict:

    command['modules'].append('assemble')
    command['params'].append(f"--assembler {kwargs['assembler']}")
    if kwargs["spades_args"] != "":
        command['params'].append(f"--spades_args {kwargs['spades_args']}")
    
    return command