#!/usr/bin/env python3

import sys,pandas as pd

tool = sys.argv[1]
data = sys.argv[2]
df = pd.read_csv(data, sep="\t")

if tool == "sylph":
    sp = df["Species"].values[0]
else:
    sp = df["Match 1"].values[0]


print(sp)