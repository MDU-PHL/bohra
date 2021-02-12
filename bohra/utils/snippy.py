import toml, pathlib, subprocess, sys
from snakemake import shell



def generate_snippy_cmd(r1, r2, isolate, reference, threads):
    
    # print(r1)
    # print(r2)
    # print(isolate)
    # print(reference)
    # print(threads)
    cmd = f'snippy --outdir {isolate} --ref {reference} --R1 {r1} --R2 {r2} --force --cpus {threads}'
    print(cmd)
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
        print(f"Reads for isolate {isolate} have passed checks and snp calling will now be peformed : {cmd}.")
        p = run_cmd(cmd)
        if p == 0:
            print(f'snp calling was successful - snippy finished running and returned a 0 exit code. Phew!')
            data[isolate]['snippy']['alignment'] = f"{isolate}/snps.aligned.fa"
            data[isolate]['snippy']['vcf'] = f"{isolate}/snps.vcf"
            data[isolate]['snippy']['cmd'] = cmd
            data[isolate]['snippy']['done'] = 'Yes'
        else:
            print('Something went wrong with snp calling. It seems that snippy struggled - please check reads and reference are ok and try again.')
            data[isolate]['snippy']['done'] = 'No'
    else:
        print(f'Isolate {isolate} did not pass checks and will not be included in further analysis' )
        data[isolate]['snippy']['done'] = 'No'
    print('Saving toml file for snippy.')
    write_toml(data = data, output = f"{isolate}/snippy.toml") 

print('Hello from snippy.py :D')
inputs = snakemake.input
isolate = snakemake.wildcards.sample
output = snakemake.output
reference = snakemake.params.reference
threads = snakemake.threads

main(inputs = inputs, isolate = isolate, output = output,reference = reference, threads =threads)
