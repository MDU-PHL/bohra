import toml, pathlib, subprocess, sys, pandas
from snakemake import shell

def generate_abritamr_cmd(input_file, isolate, wd, job_id):
    
    w = pathlib.Path(wd, isolate)
    cmd = f"abriTAMR -c {input_file} -pfx {isolate} -w {w}"
    return cmd

def generate_dict(isolate):

    df = pandas.read_csv(f"{isolate}/summary_matches.csv")
    return df.to_dict(orient = 'records')

def run_cmd(cmd):
    
    p = subprocess.run(cmd, shell = True, capture_output=True, encoding = 'utf-8')
    return p.stderr

def open_toml(tml):

    data = toml.load(tml)

    return data

def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)
    
def main(inputs, isolate, seqdata, wd, job_id):
    
    seqdata = open_toml(seqdata)
    data = {}
    data[isolate] = {}
    data[isolate]['resistome'] = {}
    if seqdata[isolate]['seqdata']['data']['Quality'] == 'PASS':
    # set up data dict
        tml = open_toml(inputs)
        contigs = f"{isolate}/contigs.fa"
        cmd = generate_abritamr_cmd(input_file = contigs, isolate = isolate, wd = wd, job_id = job_id)
        p = run_cmd(cmd)
        data[isolate]['resistome']['tool'] = 'abritamr'
        
        if f"pipeline successfully completed" in p:
            data[isolate]['resistome']['data'] = generate_dict(isolate)
            data[isolate]['resistome']['done'] = True
        else:
            data[isolate]['resistome']['done'] = False
    else:
        data[isolate]['resistome']['done'] = False
        data[isolate]['resistome']['tool'] = f"AMR not performed - failed QC"
    
    write_toml(data = data, output= f'{isolate}/resistome.toml')

inputs = snakemake.input.assembly
isolate = snakemake.wildcards.sample
seqdata = snakemake.input.seqdata
wd = snakemake.params.work_dir
job_id = snakemake.params.job_id

main(inputs = inputs, isolate = isolate, seqdata = seqdata, wd = wd, job_id = job_id)