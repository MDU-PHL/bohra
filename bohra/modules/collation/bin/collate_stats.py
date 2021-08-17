#!/usr/bin/env python3
import sys, json, datetime, subprocess, pathlib
from Bio import SeqIO

# Isolate\tMatch 1\t%\tMatch 2\t%\tMatch 3\t%

STATS_TEXT = ["Isolate\tReads\tYield\tGC content\tMin len\tAvg len\tMax len\tAvg qual\tEstimated average depth\tQuality"]


def get_vals_seqtk(f):

    with open(f, 'r') as inputfile:
        lines = inputfile.read().split('\n')
        qscore = lines[1].split('\t')[3]
        bases = lines[1].split('\t')[0]
        gc= float(lines[1].split('\t')[1]) + float(lines[1].split('\t')[2])
        return float(qscore), int(bases), round(gc,1)

def sum_stats(f,cmd,col):

    p = subprocess.run(f"cat {f} | datamash {cmd} --header-in {col}", shell = True, capture_output = True, encoding= 'utf-8')

    return p.stdout

def get_vals_kit(f):

    reads = sum_stats(f, 'sum', 4)
    minval = sum_stats(f, 'mean', 6)
    avgval = sum_stats(f, 'mean', 7)
    maxval = sum_stats(f, 'mean', 8)

    return int(reads.strip()), float(minval.strip()), float(avgval.strip()), float(maxval.strip())

def get_length(ref):

    p = subprocess.run(f"any2fasta {ref} | seqkit stats -a -T | cut -f5 | sed 1d", shell = True, capture_output = True, encoding = "utf-8")
    length = int(p.stdout.strip())
    print(p)
    return length

def get_dpth(ref,bases):
    length = get_length(ref)
    dpth = int(bases)/length
    
    return round(dpth, 1)

def get_qual(min_qscore, min_dpth, qscore, dpth):

    if dpth > min_dpth and qscore > min_qscore:
        return 'PASS'
    else:
        return 'FAIL'


qscore, bases, gc = get_vals_seqtk(f = sys.argv[2])
reads, minval, avgval, maxval = get_vals_kit(f = sys.argv[3])
dpth = get_dpth(ref = sys.argv[4], bases = bases)
qual = get_qual(min_qscore=float(sys.argv[5]), qscore=float(qscore), min_dpth=float(sys.argv[6]), dpth = dpth)
ROW = [sys.argv[1], reads, bases, gc, minval, avgval, maxval, qscore, dpth, qual]
ROW = [f"{r}" for r in ROW]
ROW = '\t'.join(ROW)
STATS_TEXT.append(ROW)
output = pathlib.Path(sys.argv[7])
output.write_text('\n'.join(STATS_TEXT))
# print(ROW)




    


