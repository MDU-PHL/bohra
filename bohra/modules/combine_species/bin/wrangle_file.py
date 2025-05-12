#!/usr/bin/env python3

import sys,pandas as pd

_files = {
    "reads": sys.argv[1],
    "asm": sys.argv[2]
}

reads_result = sys.argv[1]
asm_result = sys.argv[2]

def add_input_name(df, name):
    for col in df.columns:
        if col != "Isolate":
            df = df.rename(columns = {f"{col}": f"{col} ({name})"})
        
    return df

tab = pd.DataFrame()

for _file in _files:
    if _files[_file] != "no_results":
        df = pd.read_csv(_files[_file], sep="\t")
        df = add_input_name(df, _file)
        if tab.empty:
            tab = df
        else:
            tab = pd.merge(tab, df, on="Isolate", how="outer")
tab.to_csv(sys.stdout, sep="\t", index=False)