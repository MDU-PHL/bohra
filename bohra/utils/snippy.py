import toml, pathlib, subprocess, sys
from snakemake import shell

def generate_snippy_cmd(r1, r2, isolate, reference, threads):
    
    p = pathlib.Path(isolate)
    cmd = f"snippy --outdir {isolate} --ref {reference} --R1 {r1} --R2 {r1} --force --cpus {threads}"
    
    return cmd

def run_cmd(cmd):
    
    p = subprocess.run(cmd, shell = True, capture_output=True, encoding = 'utf-8')

    return p.returncode

def get_reads(inputs, isolate):

    s = open_toml(inputs)
    r1 = s[isolate]['seqdata']['R1']
    r2 = s[isolate]['seqdata']['R2']
    return r1,r2

def get_quality(inputs, isolate):

    s = open_toml(inputs)
    if s[isolate]['seqdata']['data']['Quality'] == 'PASS':
        return 'Yes'
    else:
        return 'No'


def open_toml(tml):

    data = toml.load(tml)

    return data

def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)
    
def main(inputs, isolate, output, reference, threads):
    
    r1,r2 = get_reads(inputs = inputs, isolate = isolate)
    run_snippy = get_quality(inputs = inputs, isolate = isolate)
    data = {}
    data[isolate] = {}
    data[isolate]['snippy'] = {}
    data[isolate]['snippy']['reference'] = reference
    data[isolate]['snippy']['R1'] = r1
    data[isolate]['snippy']['R2'] = r2
    data[isolate]['snippy']['run_snippy'] = run_snippy

    if run_snippy == 'Yes':
        cmd = generate_snippy_cmd(r1 = r1, r2=r2, isolate = isolate, reference = reference, threads = threads)
        p = run_cmd(cmd)
        if p == 0:
            data[isolate]['snippy']['alignment'] = f"{isolate}/snps.aligned.fa"
            data[isolate]['snippy']['vcf'] = f"{isolate}/snps.vcf"
            data[isolate]['snippy']['cmd'] = cmd
            data[isolate]['snippy']['done'] = 'Yes'
        else:
            data[isolate]['snippy']['done'] = 'No'
    write_toml(data = data, output = f"{isolate}/snippy.toml") 

inputs = snakemake.input
isolate = snakemake.wildcards.sample
output = snakemake.output
reference = snakemake.params.reference
threads = snakemake.threads

main(inputs = inputs, isolate = isolate, output = output,reference = reference, threads =threads)
