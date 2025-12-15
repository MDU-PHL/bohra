#!/usr/bin/env python3


import pandas as pd
import sys
import numpy as np
import argparse


def _assess(args):
    
    qc = []
    for f in args.snippy_qc:
        try:
            df = pd.read_csv(f, sep='\t')
            qc.append(df)
        except Exception as e:
            # print(f"Error reading {f}: {e}")
            continue

    qc = pd.concat(qc, ignore_index=True)
    sd2 = qc["% Aligned"].std()*2
    mean = qc["% Aligned"].mean()
    lower = mean - sd2
    
    qc["Aln_outlier"] = np.where(qc["% Aligned"] < lower, f"Alignment percentage is below the expected range of {lower}.", "")
    qc.to_csv(args.output, sep='\t', index=False)
    print(lower)


def set_parsers():
    # setup the parser
    parser = argparse.ArgumentParser(description='Collate and write bohra report',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--snippy_qc',
        help='',
        nargs='+',
        default = '')
    parser.add_argument('--output',
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





