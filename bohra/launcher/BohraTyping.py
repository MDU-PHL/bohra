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

def _assembly_required(input_table: pd.DataFrame, command:dict, kwargs:dict) -> dict:
    """Check if assembly is required based on the provided arguments."""
    LOGGER.info("Checking if assembly is required.")
    
    ctg_list = [i for i in input_table['assembly'].unique().tolist() if i != 'no_contigs']

    if 'no_contigs' in input_table['assembly'].unique().tolist() and len(ctg_list) > 0:
        LOGGER.info("Assembly is required in order to run the pipeline")
        command = _setup_assembly_args(kwargs, command)
    else:
        LOGGER.info("Assembly is not required.")
    
    if _missing_reads(input_table):
        LOGGER.warning("Some reads are missing from the input file, were assemblies are not present. This may be an error. Note that the typing pipeline requires assemblies, you may get unexpected results.")
        
        
    return command
        
def _setup_typing_args(kwargs:dict, command:dict) -> dict:

    input_table = _open_input_file(kwargs['input_file'])
    command = _assembly_required(input_table = input_table, command = command, kwargs=kwargs)
    command['modules'].append('typing')
    
    return command