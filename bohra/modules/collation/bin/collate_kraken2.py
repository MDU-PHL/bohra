#!/usr/bin/env python3
import sys, json, datetime, subprocess, pathlib
from subprocess import PIPE
import pandas as pd
import numpy as np

# Isolate\tMatch 1\t%\tMatch 2\t%\tMatch 3\t%

# SPECIES_TEXT = ["Isolate\tMatch 1\tDetected_1 (%)\tMatch 2\tDetected_2 (%)\tMatch 3\tDetected_3 (%)"]
# RESULT = [f"{sys.argv[1]}"]
# p = subprocess.run(f"grep -P '\tS\t' {sys.argv[2]}  | head -n 3" , shell = True, stdout=PIPE, stderr=PIPE, encoding = 'utf-8')

# species = p.stdout.split('\n')
# for s in species:
#     ln = s.split('\t')
#     if len(ln) == 6:
#         RESULT.append(ln[5].strip())
#         RESULT.append(ln[0].strip())

# SPECIES_TEXT.append("\t".join(RESULT))

# outfile = pathlib.Path(sys.argv[3])

# outfile.write_text('\n'.join(SPECIES_TEXT))

try:
    sp = pd.read_csv(f"{sys.argv[3]}", sep = "\t", header = None, names = ["perc","cumulative",'direct','rank','taxon_id','name'])
    sp['name'] = sp['name'].apply(lambda x: x.strip())
    uc = sp[sp['rank'] == "U"]
    host= sp[sp['name'] == sys.argv[2]]
    if uc.shape[0] == 0:
        uc = 0
    else:
        uc = uc['perc'].values[0]
    if host.shape[0] == 0:
        hst = 0
    else:
        hst = host['perc'].values[0]
    print(uc)
    spc = sp[sp['rank'] == "S"].sort_values(['perc'], ascending = False).head(3)
    print(spc)

    levs = 1
    res = {
        "Isolate": f"{sys.argv[1]}",
        "Unclassified (%)" : uc,
        f"{sys.argv[2]} (%)": hst
        }
    for row in spc.iterrows():
        res[f"Match {levs}"] = row[1]["name"]
        res[f"Detected_{levs} (%)"] = row[1]["perc"]
        levs = levs +1
    
    df = pd.DataFrame(res, index = [0])
    cols = ["Isolate","Unclassified (%)",f"{sys.argv[2]} (%)", "Match 1","Detected_1 (%)","Match 2","Detected_2 (%)","Match 3","Detected_3 (%)"]
    cols = [i for i in cols if i in df.columns.tolist()]
    print(df[cols])
    df[cols].to_csv(f"{sys.argv[4]}", sep = "\t", index = False)
except Exception as e:
    print(f"Something has gone wrong opening the kraken report: {e}")

    