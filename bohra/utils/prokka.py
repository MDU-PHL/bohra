import toml, pathlib, subprocess, sys
from snakemake import shell

def generate_prokka_cmd(isolate, assembly, threads):
    
    cmd =  f"prokka --outdir {isolate} --prefix {isolate} --mincontiglen 500 --notrna --fast --force {assembly} --cpus {threads}"

    return cmd

def generate_rm_cmd(isolate):

    cmd = f"rm {isolate}/*.err {isolate}/*.faa {isolate}/*.ffn {isolate}/*.fsa {isolate}/*.sqn {isolate}/*.tbl {isolate}/*.tsv"
    return cmd

def run_cmd(cmd):
    
    print(f'Running : {cmd}')
    p = subprocess.run(cmd, shell = True, capture_output=True, encoding = 'utf-8')
    return p.returncode
   

def open_toml(tml):

    data = toml.load(tml)

    return data

def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)
    
def main(inputs, isolate, seqdata, threads):
    
    # set up data dict
    data = open_toml(inputs)
    # data[isolate] = {}
    data[isolate]['prokka'] = {}
    # run kraken
    assembly = f"{isolate}/contigs.fa"
    seqdata = open_toml(seqdata)
    cmd = generate_prokka_cmd(isolate = isolate, assembly = assembly, threads = threads)
    if seqdata[isolate]['seqdata']['data']['Quality'] == 'PASS':
        print(f"Isolate {isolate} has passed quality checks, prokka will be used to annotate assembly.")
        p = run_cmd(cmd)
        if p == 0:
            print(f"Prokka ran successfully, now cleaning up.")
            rm_cmd = generate_rm_cmd(isolate = isolate)
            r = generate_rm_cmd(run_cmd)
            data[isolate]['prokka']['done'] = True
            data[isolate]['prokka']['gff'] = f'{isolate}/{isolate}.gff'
            data[isolate]['prokka']['txt'] = f'{isolate}/{isolate}.txt'
    else:
        print(f"Isolate {isolate} did not pass quality checks no further analysis will be performed.")
        data[isolate]['prokka']['done'] = False

    
    write_toml(data = data, output= f'{isolate}/prokka.toml')

inputs = snakemake.input.assembly
isolate = snakemake.wildcards.sample
seqdata = snakemake.input.seqdata
threads = snakemake.threads
main(inputs = inputs, isolate = isolate, seqdata = seqdata, threads = threads)
    

