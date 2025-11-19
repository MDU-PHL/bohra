from bohra.launcher.Utils import CustomFormatter, _check_path, _check_size_file
from bohra.launcher.BohraSetupFiles import _open_input_file, _check_sequence_file
from bohra.launcher.BohraAssembly import  _setup_assembly_args
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


def _check_databases(path :str, dtbtype:str) -> dict:
    
    if path == "":
        LOGGER.info(f"No database path provided using default databases installed for {dtbtype}.")
        
        return f"--{dtbtype} no_db"
    elif _check_path(path):
        return f"--{dtbtype} {path}"


def _missing_reads(input_table: pd.DataFrame) -> bool:

    noctg = input_table[input_table['assembly'] == 'no_contigs']
    if noctg.empty:
        LOGGER.info("There is nothing missing from the input file.")
        res = False
    elif not noctg.empty:
        res = False
        for row in noctg.iterrows():
            if not _check_sequence_file(pathlib.Path(row[1]['r1'])):
                LOGGER.warning(f"Missing reads for {row[1]['Isolate']}.")
                res =  True
        
    return res

def _setup_typing_args(kwargs:dict, command:dict, mtb:False) -> dict:

    typing_dbs = ["blast_db", "data_dir", "mobsuite_db"]
    input_table = _open_input_file(kwargs['input_file'])
    command = _setup_assembly_args(kwargs, command, mtb)
    command['modules'].append('typing')
    for d in typing_dbs:
        db = _check_databases(path = kwargs[d], dtbtype = d)
        command["params"].append(db)

    exclude = eval(kwargs['mlst_exclude'])
    
    if exclude != []:
        command['params'].append(f"--mlst_exclude {','.join(exclude)}")

    return command