#!/usr/bin/env python3
import sys,  json
import pandas as pd


species = sys.argv[1]

_files = sys.argv[2:]

for f in _files:
    df = pd.read_csv(f, sep='\t')
    df["Species"] = species
    # out_file = f.replace('.txt', f'_{species}_matches.txt')
    df.to_csv(f, sep='\t', index=False)