#!/usr/bin/env python3
import pathlib, subprocess, sys, pandas, numpy, csv

import pandas, pathlib


_file = sys.argv[1]

results = {'Isolate': sys.argv[2]}
with open(_file, 'r') as f:
    reader = csv.DictReader(f, delimiter = '\t')
    for row in reader:
        if row['mash_nearest_neighbor'] not in results:
            results[row['mash_nearest_neighbor']] = f"Contigs: {row['num_contigs']} - total length {row['size']}; Mash_distance: {row['mash_neighbor_distance']}"
            

header = []

df = pandas.DataFrame(results, index = [0])
df.to_csv('plasmid.txt', sep = '\t', index = False)
