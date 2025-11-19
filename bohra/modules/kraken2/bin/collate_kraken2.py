#!/usr/bin/env python3
import sys, json, datetime, subprocess, pathlib
from subprocess import PIPE

# Isolate\tMatch 1\t%\tMatch 2\t%\tMatch 3\t%

SPECIES_TEXT = ["Isolate\tMatch 1\tDetected_1 (%)\tMatch 2\tDetected_2 (%)\tMatch 3\tDetected_3 (%)"]
RESULT = [f"{sys.argv[1]}"]
p = subprocess.run(f"grep -P '\tS\t' {sys.argv[2]}  | head -n 3" , shell = True, stdout=PIPE, stderr=PIPE, encoding = 'utf-8')

species = p.stdout.split('\n')
for s in species:
    ln = s.split('\t')
    if len(ln) == 6:
        RESULT.append(ln[5].strip())
        RESULT.append(ln[0].strip())

SPECIES_TEXT.append("\t".join(RESULT))

outfile = pathlib.Path(sys.argv[3])

outfile.write_text('\n'.join(SPECIES_TEXT))
