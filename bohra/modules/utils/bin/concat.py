#!/usr/bin/env python3
import pathlib, pandas, math, sys,  re


def aggregate_results(result_files, output):
    result_df = []
    for r in result_files:
        try:
            tmp = pandas.read_csv(r, sep = '\t')
            result_df.append(tmp)
        except:
            print(f"{r} has no columns")
    if result_df == []:
        print("No valid result files found")
        sys.exit(1)
    result_df = pandas.concat(result_df)
    result_df = result_df.fillna('')

    result_df.to_csv(f"{output}.txt", sep = '\t', index = False)


result_files = []

for i in sys.argv[2:]:
    l = i.strip('[,]')
    # print(l)
    result_files.append(l)

aggregate_results(result_files = result_files, output = sys.argv[1])