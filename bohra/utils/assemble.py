import toml, pathlib, subprocess, sys, os

from snakemake import shell

def generate_asm_cmd(assembler, r1, r2, isolate, threads = 4, memory = 8):
    
    if assembler == 'shovill':
        cmd = f"shovill --R1 {r1} --R2 {r2} --outdir {isolate}/shovill --cpus {threads} --ram {memory} --force"
    elif assembler == 'spades':
        cmd = f"spades.py -1 {r1} -2 {r2} -o {isolate}/spades --threads {threads} --memory {memory}"
    else:
        cmd = f"skesa --fastq {r1},{r2} --contigs_out {isolate}/contigs.fa --cores {threads} --memory {memory}"
    print(f"Assembling using : {cmd}")
    return cmd


def generate_cmd(prefill, r1, r2, isolate, assembler, data):
    
    prfl = pathlib.Path(prefill, isolate, 'contigs.fa')s
    if prfl.exists() and os.access(f"{prfl}", os.R_OK):
        cmd = f"cp {prfl} {isolate}/contigs.fa"
        data[isolate]['assembly']['source'] = f"Prefilled from: {prfl}"
        print(f'Isolate has passed quality checks and assembly will be retrieved from {prfl}.')
    else:
        print(f'Isolate has passed quality checks and assembly will be generated.')
        cmd = generate_asm_cmd(assembler = assembler, r1 = r1, r2 = r2, isolate = isolate)
        data[isolate]['assembly']['source'] = cmd
    
    return cmd, data

def generate_mv_cmd(assembler, isolate):

    cmd = f"mv {isolate}/{assembler}/contigs.fa {isolate}/contigs.fa"
    return cmd

def run_cmd(cmd):
    
    p = subprocess.run(cmd, shell = True, capture_output=True, encoding = 'utf-8')
    return p.returncode
   

def open_toml(tml):

    data = toml.load(tml)

    return data

def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)
    
def main(seqdata, isolate, assembler, prefill):
    
    # set up data dict
    seqdata = open_toml(tml = seqdata)
    r1 = seqdata[isolate]['seqdata']['R1']
    r2 = seqdata[isolate]['seqdata']['R2']
    data = {}
    data[isolate] = {}
    data[isolate]['assembly'] = {}
    # run kraken
    # data[isolate]['seqdata']['data']['Quality']
    print(f'Checking the quality of isolate {isolate}')
    if seqdata[isolate]['seqdata']['data']['Quality'] == 'PASS':
        
        cmd, data = generate_cmd(prefill = prefill, r1 = r1, r2 = r2, isolate = isolate, assembler =assembler, data = data)
        assembly_returncode = run_cmd(cmd)
        if assembly_returncode == 0:
            print('Assembly was successful, toml will be updated.')
            p = run_cmd(generate_mv_cmd(assembler = assembler, isolate = isolate))
            # add to data dict
            data[isolate]['assembly']['done'] = 'Yes'
            data[isolate]['assembly']['assembler'] = assembler
    else:
        print(f'Isolate {isolate} did not pass quality checks so assembly was not performed.')
        data[isolate]['assembly']['done'] = 'No'
        data[isolate]['assembly']['assembler'] = f"Assembly not performed - failed QC"

    write_toml(data = data, output= f'{isolate}/assembly.toml')

seqdata = snakemake.input
isolate = snakemake.wildcards.sample
assembler = snakemake.params.assembler
prefill= snakemake.params.prefill_path

main(seqdata = seqdata, isolate = isolate,  assembler = assembler, prefill = prefill)