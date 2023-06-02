#!/usr/bin/env python3
import json,sys,pathlib, subprocess


p = subprocess.run(f"grep -P '\tS\t' {sys.argv[1]}  | head -n 1" , shell = True, capture_output = True, encoding = 'utf-8')

species = p.stdout.split('\n')[0].split('\t')[-1].strip()
print(species)