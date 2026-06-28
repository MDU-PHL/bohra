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
fh = logging.FileHandler('bohra.log')
fh.setLevel(logging.DEBUG)
formatter = logging.Formatter('[%(levelname)s:%(asctime)s] %(message)s', datefmt='%Y-%m-%d %I:%M:%S %p') 
fh.setFormatter(formatter)
LOGGER.addHandler(ch) 
LOGGER.addHandler(fh)


def _check_databases(path :str, dtbtype:str) -> dict:
    
    if path == "":
        LOGGER.info(f"No database path provided using default databases installed for {dtbtype}.")
        return f"--{dtbtype} no_db"
    elif _check_path(path):
        LOGGER.info(f"Found database path for {dtbtype}.")
        return f"--{dtbtype} {path}"
    else:
        LOGGER.info(f"No database path provided using default databases installed for {dtbtype}.")
        return f"--{dtbtype} no_db"

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

    # typing_dbs = ["mlst_dbdir", "mobsuite_db"]
    # LOGGER.info(f"adding dbs")
    typing_dbs = {
        "blast_db": f"{kwargs['mlst_dbdir']}/blast/mlst.fa",
        "data_dir":f"{kwargs['mlst_dbdir']}/pubmlst",
        "mobsuite_db": f"{kwargs['mobsuite_db']}"}
    input_table = _open_input_file(kwargs['input_file'])
    command = _setup_assembly_args(kwargs, command, mtb)
    command['modules'].append('typing')
    for d in typing_dbs:
        LOGGER.info(f"Checking the {d} db at {typing_dbs[d]}")
        db = _check_databases(path = typing_dbs[d], dtbtype = d)
        # LOGGER.info(f"{db} will be added to run paramaters")
        command["params"].append(db)

    exclude = eval(kwargs['mlst_exclude'])
    
    if exclude != []:
        command['params'].append(f"--mlst_exclude {','.join(exclude)}")

    return command