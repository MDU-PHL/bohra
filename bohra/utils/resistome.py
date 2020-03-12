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
    
    print(f"Running : {cmd}.")
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
        print(f"Isolate {isolate} passed quality checks, resistome will be found.")
        tml = open_toml(inputs)
        contigs = f"{isolate}/contigs.fa"
        cmd = generate_abritamr_cmd(input_file = contigs, isolate = isolate, wd = wd, job_id = job_id)
        p = run_cmd(cmd)
        data[isolate]['resistome']['tool'] = 'abritamr'
        
        if f"pipeline successfully completed" in p:
            print(f"Resistome was successfully identified, results will now be added to the toml file.")
            data[isolate]['resistome']['data'] = generate_dict(isolate)
            data[isolate]['resistome']['done'] = True
        else:
            print(f"Something has gone wrong with the resistome - please check your data and you may try again.")
            data[isolate]['resistome']['done'] = False
    else:
        print(f"Isolate {isolate} did not pass quality checks - no further analysis will be performed.")
        data[isolate]['resistome']['done'] = False
        data[isolate]['resistome']['tool'] = f"AMR not performed - failed QC"
    
    write_toml(data = data, output= f'{isolate}/resistome.toml')

inputs = snakemake.input.assembly
isolate = snakemake.wildcards.sample
seqdata = snakemake.input.seqdata
wd = snakemake.params.work_dir
job_id = snakemake.params.job_id

main(inputs = inputs, isolate = isolate, seqdata = seqdata, wd = wd, job_id = job_id)