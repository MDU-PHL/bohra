#!/usr/bin/env python3
import sys, json, datetime, subprocess, pathlib, pandas

files = sys.argv[2:]
df = pandas.DataFrame()
for _file in files:
    print(_file)
    tmp = pandas.read_csv(_file, sep = '\t')
    if df.empty:
        df = tmp
    else:
        df = df.append(tmp)
df = df.fillna('')
df.to_csv(f'{sys.argv[1]}.txt', sep = '\t', index = False)
    
