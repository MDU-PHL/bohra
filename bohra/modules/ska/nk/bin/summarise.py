#!/usr/bin/env python3

import pandas as pd
import sys



pth = sys.argv[1]
out_stats = sys.argv[2]

with open(pth, 'r') as f:
    lines = f.read().strip().split('\n')

union = int([line for line in lines if line.startswith('k-mers=')][0].split('=')[1].strip())
names = eval([line for line in lines if line.startswith('sample_names=')][0].split('=')[1])
numkmers = eval([line for line in lines if line.startswith('sample_kmers=')][0].split('=')[1])

raw = dict(zip(names, numkmers))

data = []

for r in raw:
    d = {
        'Isolate': r,
        '# kmers': raw[r],
        'Union # kmers': union,
        'Proportion in union (%)': round((raw[r] / union) * 100, 2)
    }

    data.append(d)

df = pd.DataFrame(data)
df.to_csv(out_stats, sep = '\t', index = False)