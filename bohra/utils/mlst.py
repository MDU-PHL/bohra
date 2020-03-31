import toml, pathlib, subprocess, sys, pandas, json
from snakemake import shell

def generate_mlst_cmd(assembly, isolate):

    cmd = f"mlst --label {isolate} --json {isolate}/{isolate}.json {assembly}"
    return cmd

def extract_mlst(isolate):

    jsn = json.load(open(f"{isolate}/{isolate}.json"))
    if isinstance(jsn, list):
        pos = 0
    else:
        pos = isolate
    # isolate = jsn[0]['id']
    alleles = []
    if jsn[pos]['alleles']:
        ax = zip(jsn[pos]['alleles'].keys(), jsn[pos]['alleles'].values())
        for a in ax:
            alleles.append(f"{a[0]}({a[1]})")
        alleles = sorted(alleles)

    d = {
        'Isolate':isolate,
        'alleles': alleles,
        'len_alleles': len(alleles),
        'scheme':jsn[pos]['scheme'],
        'ST':jsn[pos]['sequence_type']
    }
    return d

def generate_rm_cmd(isolate):

    cmd = f"rm {isolate}/{isolate.json}"
    return df.to_dict(orient = 'records')

def run_cmd(cmd):
    print(f"Running : {cmd}")
    p = subprocess.run(cmd, shell = True, capture_output=True, encoding = 'utf-8')
    return p.returncode

def open_toml(tml):

    data = toml.load(tml)

    return data

def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)
    
def main(inputs, isolate, seqdata):
    
    print(f"Hello from mlst. Inputs are : {inputs}")
    data = {}
    data[isolate] = {}
    data[isolate]['mlst'] = {}
    seqdata = open_toml(seqdata)
    assembly = f"{pathlib.Path(isolate, 'contigs.fa')}"
    if seqdata[isolate]['seqdata']['data']['Quality'] == 'PASS':
    # set up data dict
        print(f"Isolate {isolate} has passed quality checks, mlst will be found.")
        cmd = generate_mlst_cmd(assembly = assembly, isolate = isolate)
        p = run_cmd(cmd)
        if p == 0:
            print('MLST has successfuly completed, toml will now be updated.')
            data[isolate]['mlst']['done'] = 'Yes'
            data[isolate]['mlst']['data'] = extract_mlst(isolate)
    else:
        print(f"Isolate {isolate} did not pass quality checks, no mlst will be determined.")
        data[isolate]['mlst']['done'] = 'No'
        data[isolate]['mlst']['data'] = f"MLST not performed - failed QC"
        
    write_toml(data = data, output= f'{isolate}/mlst.toml')

inputs =snakemake.input.assembly
isolate = snakemake.wildcards.sample
seqdata = snakemake.input.seqdata

main(inputs = inputs, isolate = isolate, seqdata = seqdata)
    

