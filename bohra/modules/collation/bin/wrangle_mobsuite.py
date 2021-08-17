#!/usr/bin/env python3
import pathlib, subprocess, sys, pandas, numpy, csv

import pandas, pathlib


_file = sys.argv[1]

results = {'Isolate': sys.argv[2]}
with open(_file, 'r') as f:
    reader = csv.DictReader(f, delimiter = '\t')
    for row in reader:
        if row['molecule_type'] == 'plasmid':
            if row['mash_nearest_neighbor'] not in results:
                results[row['mash_nearest_neighbor']] = f"{row['contig_id']}"
            else:
                results[row['mash_nearest_neighbor']] = f"{results[row['mash_nearest_neighbor']]},{row['contig_id']}"

header = []

df = pandas.DataFrame(results, index = [0])
df.to_csv('plasmid.txt', sep = '\t', index = False)
