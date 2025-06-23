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
            return False
    else:
        return False





def _setup_basic_args(kwargs:dict, command:dict) -> dict:

    mods = command.get('modules', [])
    params = command.get('params', [])
    spn = _is_speciation(kwargs["speciation"])
    if spn:
        mods.append(spn)
        params.append(_species_tool(kwargs["speciation"]))
        chk_db = _check_species_database(kwargs["speciation"], kwargs[f"{kwargs['speciation']}_db"])
        if chk_db:
            params.append(chk_db)
    
    
    command['modules'] = mods
    command['params'] = params

    return command