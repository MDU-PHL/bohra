#!/usr/bin/env python3
import sys, json, datetime, subprocess, pathlib,pandas

# Isolate\tMatch 1\t%\tMatch 2\t%\tMatch 3\t%

STATS_TEXT = ["Isolate\tReads\tYield\tGC content\tMin len\tAvg len\tMax len\tAvg qual\tEstimated average depth\tQuality (>Q30)"]


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

# STATS_TEXT = ["Isolate\tReads\tYield\tGC content\tMin len\tAvg len\tMax len\tAvg qual\tEstimated average depth\tQuality (>Q30)"]
gcs = pandas.read_csv(sys.argv[3], sep = '\t')
tab = pandas.read_csv(sys.argv[2], sep = '\t')
tab['Isolate'] = sys.argv[1]
print(tab)
print(gcs)
# file   format   type   num_seqs   sum_len     min_len   avg_len   max_len   Q1      Q2      Q3      sum_gap   N50   Q20(%)   Q30(%)
tab = tab.rename(columns = {'num_seqs': 'Reads', 'sum_len': 'Yield','min_len':'Min len','max_len':'Max len', 'avg_len':'Avg len', 'Q30(%)': 'Average quality (% >Q30)'})
tab['Estimated average depth'] = get_dpth(ref = sys.argv[4], bases = tab['Yield'].values[0])
tab['GC content'] = gcs[gcs.columns[0]].values[0]
tab['Avgerage quality'] = gcs[gcs.columns[1]].values[0]
tab = tab[['Isolate','Reads','Yield','GC content','Min len','Avg len','Max len','Average quality (% >Q30)','Estimated average depth']]
tab.to_csv('read_assessment.txt', sep = '\t', index = False)
# 
# print(ROW)




    


