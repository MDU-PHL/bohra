#!/usr/bin/env python3
import pathlib, subprocess, sys, pandas

from Bio import SeqIO
                
HEADER = ["Isolate\tLength\tAligned\tUnaligned\tHeterozygous\tLow coverage\t% Aligned"]

def check_snippy(isolate,aln):
    
    p = pathlib.Path(aln)
    fasta = p.open()
    length = 0
    nocov = 0
    lowcov = 0
    het = 0
    unaln = 0
    for i in SeqIO.parse(fasta,'fasta'): # use BioPython to determine percent alignment
        length = length + len(i.seq)
        nocov = nocov + i.seq.count('-')
        lowcov = lowcov + i.seq.count('N')
        het = het + i.seq.count('n')
    unaln = unaln + nocov + lowcov + het
    aln = length - unaln
    perc_aln = 100*(length - unaln) / length
    # quality = 'PASS' if perc_aln > float(minaln) else 'FAIL'
    return length,aln, nocov, lowcov, het, unaln, perc_aln

def main(inputs, isolate, output):

    length, aln, nocov, lowcov, het, unaln, perc_aln= check_snippy(aln = inputs, isolate = isolate)
    data = [isolate, length, aln, unaln,het,lowcov,round(perc_aln, 2)]
    data = [f"{d}" for d in data]
    data = "\t".join(data)
    HEADER.append(data)
    out = pathlib.Path(output)
    out.write_text("\n".join(HEADER))
    # print(data)
   
# {input} {wildcards.sample} {output} {params.minaln}
inputs = sys.argv[2]
isolate = sys.argv[1]
output = sys.argv[3]
minaln= sys.argv[4]

main(inputs = inputs, isolate = isolate, output = output)
