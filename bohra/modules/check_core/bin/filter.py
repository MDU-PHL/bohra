#!/usr/bin/env python3


import pandas as pd
import sys
import numpy as np
import argparse


def _assess(args):
    
    qc = pd.read_csv(f"{args.qc}", sep = "\t")
    qc = qc.fillna("")
    if args.strict == 'true':
        fail = [i for i in qc["Aln_outlier"].unique() if "below" in i]

        if fail != []:
            print("exclude")
        else:
            print("include")
    else:
        qc = qc[(qc["Isolate"] == args.seqid) & (qc['Aln_outlier'].str.contains("below"))]
        if qc.empty:
            print("inlcude")
        else:
            print("exclude")



def set_parsers():
    # setup the parser
    parser = argparse.ArgumentParser(description='Collate and write bohra report',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--qc',
        help='',
        default = '')
    parser.add_argument('--seqid',
        help='',
        default = '')
    parser.add_argument('--strict',
        help='',
        default = '')
    
    parser.set_defaults(func=_assess)
    args = parser.parse_args()
    
    if vars(args) == {}:
        parser.print_help(sys.stderr)
    else:
        # print(args)
        args.func(args)
	

def main():
    """
    run pipeline
    """

    args = set_parsers()
    

if __name__ == "__main__":
    main()





