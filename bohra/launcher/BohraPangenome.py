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






def _setup_pangenome_args(kwargs:dict, command:dict, mtb:False) -> dict:

    command['modules'].append('pangenome')
    command["params"].append(f"--pangenome_groups {kwargs['pangenome_groups']}")
    return command