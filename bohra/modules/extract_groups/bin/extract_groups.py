#!/usr/bin/env python3

import sys, json, pathlib
import pandas as pd





def extract_unclustered(df, col, col1):

    u = df.groupby(col)[col1].agg('count').reset_index()

    uc = list(u[u[col1] == 1][col])

    return uc

groups = sys.argv[1]



df = pd = pd.read_csv(groups, sep="\t")

if "ST" in df.columns:
    cols = [i for i in df.columns if "Allele" in i]
    df["MLST_alleles"] = df[cols].apply(lambda x: ":".join(x), axis=1)

    uc = extract_unclustered(df, col="MLST_alleles", col1="Isolate")

    df = df[~df["Isolate"].isin(uc)]

    df[["ID","MLST_alleles"]].to_csv("group.txt", sep="\t", index=False, header = False)
else:
    cols = sorted([i.split(":")[-1] for i in df.columns if "Tx" in i])

    df = df[df[f"Tx:{cols[-1]}"] != "UC"]

    df[["ID", f"Tx:{cols[-1]}"]].to_csv("group.txt", sep="\t", index=False, header = False)