#!/usr/bin/env python3
import pathlib, pandas, math, sys,  re, argparse,json,os,time
def get_date_sequenced(reads):

    ti_c = os.path.getctime(reads)
    ti_c = time.ctime(ti_c)
    t_obj = time.strptime(ti_c)

    seqd = time.strftime("%Y-%m-%d", t_obj)
    return seqd


print(get_date_sequenced(reads = sys.argv[1]))