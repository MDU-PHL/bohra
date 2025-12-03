from bohra.launcher.Utils import CustomFormatter, _check_path, _run_subprocess, _get_required_columns, _get_columns_list

import pandas as pd
import pathlib
import os
import logging
import shutil

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



def _open_input_file(_file:str) -> pd.DataFrame:
        
    """
    Open input file and return a pandas dataframe
    :param _file: path to input file
    :return: pandas dataframe
    """
    LOGGER.info(f"Opening input file {_file}.")
    if _file == "":
        LOGGER.critical(f"Input file cannot be empty. Please try again.")
        raise SystemExit
    elif not _check_path(_file):
        LOGGER.critical(f"Input file {_file} does not exist or is not accessible.")
        raise SystemExit
    else:
        LOGGER.info(f"Backing up original file {_file} to {pathlib.Path(_file)}.bak.")
        shutil.copy(_file, f"{pathlib.Path(_file)}.bak")
        return pd.read_csv(_file, engine='python', sep = None, dtype = str)

def _check_data_format(df:pd.DataFrame) -> pd.DataFrame:
    """
    Check if the input file is in the correct format
    :param df: pandas dataframe
    :return: pd.DataFrame
    """
    
    columns_req = _get_required_columns()
    columns = _get_columns_list()
    # check if all required columns are present
    for col in columns_req["required"]:
        if col not in df.columns.tolist():
            LOGGER.critical(f"Column {col} is missing from input file.")
            raise SystemExit
    # check if all must have columns are present
    mchk = 0
    for col in columns_req["must_have"]:
        if col in df.columns.tolist():
            mchk = mchk + 1
    if mchk == 0:
        LOGGER.critical(f"Must have at least one of {','.join(columns_req['must_have'])} in your input file.")
        raise SystemExit
    # check if all optional columns are present
    for col in columns:
        if col not in df.columns:
            df[col] = "not_supplied"
    
    for col in df.columns:
        if col not in columns:
            columns.append(col)
    df = df[columns]
    df = df.fillna("not_supplied")
    LOGGER.info(f"Saving input file with columns {','.join(columns)} to {pathlib.Path('input_checked.tsv')}.")
    df.to_csv('input_checked.tsv', sep='\t', index=False, header=True)

    return df

def _check_sequence_file(_file:str) -> bool:
    """ Check if the sequence file exists and is accessible 
    :param _file: path to sequence file
    :return: boolean
    """
        
    if _file.exists() and os.access(_file, os.R_OK):
        
        LOGGER.info(f"Found {_file.name}.")
        return True
    else:
        LOGGER.warning(f"The {_file.name} does not exist or is not accessible.")
        return False

def _get_input_types(cols_provided: list, cols_required:list) -> dict:

    """
    Get the input types from the columns provided
    :param cols_provided: list of columns provided in the input file
    :return: list of input types
    """
    input_types = {}
    for c in cols_required:
        if c == "assembly":
            input_types[c] = "contigs.fa"
        elif c == "r1":
            input_types[c] = "R1.fastq.gz"
        elif c == "r2":
            input_types[c] = "R2.fastq.gz"
    return input_types

def _infer_correct_read(file_path:str) -> str:

    if '_R1' in file_path or '_1.fastq' in file_path or 'r1' in file_path:
        return "R1.fastq.gz"
    elif '_R2' in file_path or '_2.fastq' in file_path or 'r2' in file_path:
        return "R2.fastq.gz"
    else:
        return "unknown"

def _check_outdir(workdir:str, report_outdir:str, replace_report:bool = False) -> bool:

    out = pathlib.Path(workdir,report_outdir)
    if out.exists() and not replace_report:
        LOGGER.critical(f"Output directory {out} already exists. Please choose a different output directory or use the --replace-report flag to overwrite.")
        raise SystemExit
    elif out.exists() and replace_report:
        LOGGER.warning(f"Output directory {out} already exists. Overwriting as per user request.")
        shutil.rmtree(out)
        out.mkdir(parents=True, exist_ok=True)
        return True
    return True

def _make_workdir(_input:pd.DataFrame, workdir:str, report_outdir:str, replace_report:bool = False) -> bool:

    columns = _get_columns_list()
    pereads = ["r1","r2"]
    if _check_path(workdir) and _check_outdir(workdir, report_outdir, replace_report):
        
        wd = pathlib.Path(workdir)
        input_table = _open_input_file(_file=_input)
        input_table = _check_data_format(df=input_table)
        
        input_types = _get_input_types(cols_provided=input_table.columns.tolist(), cols_required=columns)
        for row in input_table.iterrows():
            
            for input_type in input_types:
                user_supplied = row[1][input_type]
                
                if _check_sequence_file(pathlib.Path(user_supplied)):
                    # check that col is actually the r1 or r2
                    target_file = input_types[input_type]
                    if input_type in pereads:
                        target_file = _infer_correct_read(user_supplied)
                    if target_file == "unknown":
                        LOGGER.warning(f"Could not infer read type from {user_supplied}. Skipping link creation.")
                        continue
                    LOGGER.info(f"Linking {user_supplied} to {wd / row[1][columns[0]]}/{target_file}.")
                    try:
                        LOGGER.info(f"Creating directory {row[1][columns[0]]} in {workdir}")
                        target = pathlib.Path(wd / row[1][columns[0]] / target_file)
                        if not target.exists():
                            pathlib.Path(f"{wd / row[1][columns[0]]}").mkdir(parents=True, exist_ok=True)
                            target.symlink_to(user_supplied)
                            if target.exists():
                                LOGGER.info(f"Successfully linked {user_supplied} to {target}.")
                            else:
                                LOGGER.critical(f"{target} does not exist after linking.")
                                raise SystemExit
                        else:
                            LOGGER.warning(f"File {target} already exists. Skipping link creation.")
                    except Exception as e:
                        LOGGER.critical(f"Could not link {user_supplied} to {wd / row[1][columns[0]]}/{target_file}. Error: {e}")
                        raise SystemExit
                else:
                    LOGGER.warning(f"File {user_supplied} does not exist or is not accessible. Skipping link creation.")
    return "input_checked.tsv"
