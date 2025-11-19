#!/usr/bin/env python3

import pandas as pd
import sys


def _smoosh_results(values: list) -> str:
    """
    Smoosh the values into a single string.

    Args:
        values (list): List of values to smoosh.

    Returns:
        str: Smooshed string of values.
    """

    if values[1]== "Susceptible":
        return 'Susceptible'
    elif "resistant" in values[1].lower() :
        return f"{values[1]} ({values[0]}|{values[2]})"
    elif "not reportable" in values[1].lower():
        return f"{values[1]} ({values[0]} | not reportable)"

    return ''

def summarise_tbtamr(input_file: str, output_file: str) -> None:
    """
    Summarise the tbtamr output file and save it to a new file.

    Args:
        input_file (str): Path to the input tbtamr output file.
        output_file (str): Path to the output summary file.
    """
    df = pd.read_csv(input_file)
    df = df.rename(columns = {"Seq_id": "Isolate"})
    drs = set([ i.split('-')[0].strip() for i in df.columns if 'mech' in i])

    for dr in drs:
        df[dr] = df[[f"{dr} - mechanisms", f"{dr} - interpretation", f"{dr} - confidence"]].apply(
            lambda x: _smoosh_results(x), axis=1
        )
    
    cols = ["Isolate", "Predicted drug resistance"] + sorted(drs) + ["Phylogenetic lineage", "Db version"]

    df[cols].to_csv(output_file, sep='\t', index=False)

input_file = sys.argv[1]
output_file = sys.argv[2]

summarise_tbtamr(input_file, output_file)