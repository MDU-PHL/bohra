import toml, pathlib, subprocess, sys, pandas
from snakemake import shell

def combine_tomls(inputs, isolate):

    final_toml = {}
    final_toml[isolate] = {}
    for i in inputs:
        tml = open_toml(i)
        # print(st)
        t = list(tml[isolate].keys())[0]
        final_toml[isolate][t] = tml[isolate][t]
        # subprocess.run(f"rm {i}", shell = True, capture_output = True, encoding = "utf-8")
    
    return final_toml

def open_toml(tml):

    data = toml.load(tml)

    return data

def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)
    
def main(inputs, isolate):
    
    # set up data dict
    final_toml = combine_tomls(inputs = inputs, isolate = isolate)
    
    write_toml(data = final_toml, output= f'{isolate}/final.toml')

isolate = snakemake.wildcards.sample
inputs = snakemake.input

main(isolate = isolate, inputs = inputs)
    

