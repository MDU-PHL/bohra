from bohra.launcher.Utils import CustomFormatter, _check_path, _get_required_columns, _check_size_file, _check_input_snippy, _check_reference, _check_mask, _compartive_args
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

def _get_annotation(input_file:str, mtb:bool) -> str:
    annot_cols =[]
    if _check_path(input_file):
        df = pd.read_csv(input_file, sep='\t')
        cols = df.columns.tolist()
        std_cols = _get_required_columns()
        std = []
        for col in std_cols:
            for c in std_cols[col]:
                std.append(c)
        
        for c in cols:
            if c not in std:
                annot_cols.append(c)
    if mtb:
        annot_cols.extend(['Phylogenetic lineage', "Predicted drug resistance"])
    return ','.join(annot_cols) if annot_cols else ""

def _setup_comparative_args(kwargs:dict, command:dict, mtb:bool) -> dict:
    if mtb:
        command['modules'].append('mtb')
    command['modules'].append(kwargs['comparative_tool'])
    command['params'].append(f"--tree_input {kwargs['tree_input'] if kwargs['comparative_tool'] != 'mash' else 'distance'}")
    command['modules'].append(f"{kwargs['tree_builder'] if kwargs['comparative_tool'] != 'mash' else 'quicktree'}")
    command['params'].append(f"--cluster {'false' if kwargs['comparative_tool'] == 'mash' else kwargs['cluster']}")
    command['params'].append(f"--ignore_warnings {'true' if kwargs['ignore_warnings'] else 'false'}")
    
    if kwargs['comparative_tool'] == "snippy" and not _check_input_snippy(kwargs['input_file']):
        LOGGER.critical("You have selected snippy as your tool for snp detection. Currently this requires paired-end reads. Please check your input file and try again.")
        raise SystemExit
    # command = _accessory_params(kwargs=kwargs, command=command)
    annots = _get_annotation(kwargs['input_file'], mtb)
    ref = _check_reference(kwargs['reference_genome'])
    mask = _check_mask(kwargs['mask'])
    if annots != "":
        command['params'].append(f"--annot_cols '{annots}'")
    command['params'].append(f"--reference {ref}")
    command['params'].append(f"--mask {mask}")
    if ref == "no_ref" and kwargs['comparative_tool'] == "snippy":
        LOGGER.critical("You have selected snippy as your tool for snp detection. Currently this requires a reference genome. Please check your input file and try again.")
        raise SystemExit
    command = _compartive_args(tool = kwargs['comparative_tool'],kwargs=kwargs, command=command)
        
    if 'gubbins' in kwargs and kwargs["gubbins"]:
        command['params'].append("--gubbins true")
    
    return command