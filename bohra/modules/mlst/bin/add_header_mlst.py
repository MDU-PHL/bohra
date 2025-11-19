#!/usr/bin/env python3
import pathlib, subprocess, sys, pandas, numpy,json

import pandas, pathlib


_file = sys.argv[1]

with open(_file, 'r') as j:
    _d = json.load(j)
    header = ['Isolate','Scheme','ST']
    row = [_d[0]['id'],_d[0]['scheme'],_d[0]['sequence_type']]
    try:
        als = sorted(list(_d[0]['alleles'].keys()))
        for a in range(len(als)):
            header.append(f'Allele {a + 1}')
            row.append(f"{als[a]}({_d[0]['alleles'][als[a]]})")
    
        h = '\t'.join(header)
        r = '\t'.join(row)
    
        pathlib.Path('mlst.txt').write_text('\n'.join([h,r]))
    except:
        pathlib.Path('mlst.txt').write_text('\n'.join(['\t'.join(header), '\t'.join(row)]))