#!/usr/bin/env python3


import pandas as pd
import sys
import numpy as np
import argparse


def check_aln(aln:list)-> tuple:
    alns = [i for i in aln if not isinstance(i, str)]
    q1, q3 = np.percentile(alns, [25, 75])
    iqr = q3 - q1
    lower_bound = q1 - (1.5 * iqr)
    upper_bound = q3 + (1.5 * iqr)
    # print(f"Lower bound: {lower_bound}, Upper bound: {upper_bound}")
    return lower_bound, upper_bound



def report_aln(upper, lower, aln) -> str:
    """
    Check if the number of contigs is within the expected range
    """
    # print(contigs[0])
    # print(lower, upper)
    
    
    if (aln < lower):
        return f"Alignment percentage: {aln} is below the expected range of {lower}."
    else:
        return f""


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
    lower, upper = check_aln(qc["% Aligned"].tolist())
    qc["Aln_outlier"] = qc["% Aligned"].apply(lambda x: report_aln(upper, lower, x))
    qc.to_csv(args.output, sep='\t', index=False)
    if qc["Aln_outlier"].str.contains("below").any():
        print("outlier")
    else:
        print("ok")


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





