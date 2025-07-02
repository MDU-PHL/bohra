#!/usr/bin/env python3

import sys,pandas as pd


reads_result = sys.argv[1]
asm_result = sys.argv[2]


res = [ i for i in [reads_result.strip(), asm_result.strip()] if i != "no_results"]
# print(res)
res = list(set(res))
# print(res)
if len(res) == 0:
    print("no_result")
    sys.exit(0)
else:
    print(':'.join(res))
    sys.exit(0)
