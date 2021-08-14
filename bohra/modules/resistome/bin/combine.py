#!/usr/bin/env python3
import pathlib, pandas, math, sys,  re

def get_genes(_file, _type, _dict):
    
    with open(_file, 'r') as f:
        data = f.read().strip().split('\n')
        header = data[0].split('\t')
        line = data[1].split('\t')
        for i in range(1,len(header)):
            l = line[i] if _type != 'partials' else f"{line[i]}^"
            if header[i] in _dict:
                _dict[header[i]].append(l)
            else:
                _dict[header[i]] = [l]
        return _dict

_dict = {'Isolate': sys.argv[1]}

# matches
_dict = get_genes(_file = sys.argv[2], _type = 'matches', _dict = _dict)

# partials 
_dict = get_genes(_file = sys.argv[3], _type = 'partials', _dict = _dict)

for d in _dict:
    if d != 'Isolate':
        _dict[d] = ';'.join(_dict[d])
df = pandas.DataFrame(_dict, index = [0])
df = df.set_index('Isolate')

df.to_csv('resistome.txt', sep = '\t')
