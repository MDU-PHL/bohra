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





def _is_speciation(speciation: str) -> str:

    if speciation != "none":
        return "species"
    else:
        return False


def _species_tool(speciation:str) -> str:
    
    if speciation != 'none':
        return f"--use_{speciation} true"
    

def _check_species_database(speciation:str, database:str) -> str:
    LOGGER.info(f"Checking species database for {speciation} with path {database}")
    if _check_path(database) and speciation != "kraken2":
        return f"--{speciation}_db {database}"
    elif _check_path(database) and speciation == "kraken2":
        if pathlib.Path(database).is_dir():
                kmerfiles = sorted(pathlib.Path(database).glob('*'))
                s = []
                for k in range(len(kmerfiles)):
                    s.append(_check_size_file(pathlib.Path(database) / kmerfiles[k]))
                if 0 not in s:
                    return f"--{speciation}_db {database}"
        else:
            LOGGER.error(f"{speciation} database path: {database} is not complete. Speciation will not be performed.")
            return False
    else:
        LOGGER.error(f"{speciation} database path: {database} does not exist or is not a valid directory. Speciation will not be performed.")
        return False


def _accessory_params(kwargs:dict, command:dict) -> list:
    """Returns a list of accessory parameters for the command."""
    
    if kwargs['use_conda']:

        if _check_path(kwargs['dependency_prefix']):
            command['params'].append(f"--dependency_prefix {kwargs['dependency_prefix']}")
            command['params'].append("-with-conda")
            command['params'].append(f"--enable_conda { 'true' if kwargs['use_conda'] else 'false' }")
    if not kwargs['force']:
        command['params'].append("-resume")

    command['params'].append(f"--background_color '{kwargs['background_color']}'  --text_color '{kwargs['text_color']}' --job_id {kwargs['job_name']}")
    return command


def _setup_basic_args(kwargs:dict, command:dict) -> dict:

    
    spn = _is_speciation(kwargs["speciation"]) 
    if spn:
        chk_db = _check_species_database(kwargs["speciation"], kwargs[f"{kwargs['speciation']}_db"])
        if chk_db:
            command['params'].append(chk_db)
            command['modules'].append(spn)
            command['params'].append(_species_tool(kwargs["speciation"]))
    
    command = _accessory_params(kwargs=kwargs, command=command)

    return command