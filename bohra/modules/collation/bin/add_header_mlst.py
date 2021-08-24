#!/usr/bin/env python3
import pathlib, subprocess, sys, pandas, numpy

import pandas, pathlib


_file = sys.argv[1]

with open(_file, 'r') as f:
    alleles = len(f.read().split(',')) - 3

allele_names = ','.join([ f"Allele {i+1}" for i in list(range(alleles))])

cmd = f"csvtk add-header -T -n 'Isolate,Scheme,ST,{allele_names}' {_file} > mlst.txt"
subprocess.run(cmd, shell = True)