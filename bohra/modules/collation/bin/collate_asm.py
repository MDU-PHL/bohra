#!/usr/bin/env python3
import pathlib, subprocess, sys, pandas, numpy

import pandas, pathlib

HEADER = [f"Isolate\tAssembly length\t# Contigs\t# Gaps\tMin Contig size\tMax Contig size\tAvg Contig size\tAssembly N50\tCDS\trRNA"]

def combine(prokka, asm, isolate, output):

    gff = pandas.read_csv(prokka, sep = ':', header = None, names = ['cond', isolate])
    # print(gff)
    df = pandas.read_csv(asm, sep = '\t')
    gff = gff[gff['cond'].isin(['CDS', 'rRNA'])]
    # print(gff)
    rrna = gff[gff['cond'] == 'rRNA'][isolate].values[0].replace("\"","").strip() if not gff[gff['cond'] == 'rRNA'].empty else ''
    cds = gff[gff['cond'] == 'CDS'][isolate].values[0].replace("\"","").strip()
    # print(rrna)
    # print(cds)
    print(df)
    bp = df['sum_len'].values[0]
    contigs = df['num_seqs'].values[0]
    mincontigs = df['min_len'].values[0]
    avgcontigs = df['avg_len'].values[0]
    maxcontigs = df['max_len'].values[0]
    gaps= df['sum_gap'].values[0]
    n50 = df['N50'].values[0]

    result = f"{isolate}\t{bp}\t{contigs}\t{gaps}\t{mincontigs}\t{avgcontigs}\t{maxcontigs}\t{n50}\t{cds}\t{rrna}"
    # print(result)
    HEADER.append(result)
    out = pathlib.Path(output)
    out.write_text('\n'.join(HEADER))

def main(isolate,asm,prokka, output):
    combine(prokka=prokka, asm = asm, isolate=isolate, output = output)
    
    

isolate = sys.argv[1]
prokka = sys.argv[2]
asm = sys.argv[3]
output = sys.argv[4]
main(prokka = prokka, isolate=isolate, asm = asm, output = output)
    

