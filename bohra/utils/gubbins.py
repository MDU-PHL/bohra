import toml, pathlib, subprocess, sys, pandas, datetime
from snakemake import shell

def generate_snippy_clean_cmd():

    cmd = f"snippy-clean_full_aln core.full.aln  > clean.full.aln"								
    return cmd

def generate_gubbins_cmd():

    cmd = f"run_gubbins.py -c 36  --prefix core clean.full.aln"								
    return cmd

def generate_snp_sites_cmd():
    
    cmd = f"snp-sites -c clean.filtered_polymorphic_sites.fasta > gubbins.aln"								
    return cmd

def run_cmd(cmd):
    print(f"Now running {cmd}.")
    p = subprocess.run(cmd, shell = True, capture_output=True, encoding = 'utf-8')
    return p.returncode

def open_toml(tml):

    data = toml.load(tml)

    return data

def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)
    
def main(inputs, gubbins):
    data = {}
    data['gubbins'] = {}
    # print(gubbins)
    if gubbins != False:
        print(f"You have elected to correct for recombination. Gubbins will now be used.")
        clean = run_cmd(generate_snippy_clean_cmd())
        if clean == 0:
            gub = run_cmd(generate_gubbins_cmd())
            if gub == 0:
                sites = run_cmd(generate_snp_sites_cmd())
                if sites == 0:
                    data['gubbins']['run_gubbins'] = 'Yes'
                    data['gubbins']['aln'] = 'gubbins.aln'
    else:
        print(f"Gubbins is not necessary for the current job.")
        data['gubbins']['run_gubbins'] = 'No'

        
    write_toml(data = data, output = "gubbins.toml")

inputs = snakemake.input
gubbins = snakemake.params.gubbins
main(inputs = inputs, gubbins = gubbins)


# mash triangle -C *.msh

# mash sketch -m 5 -s 10000 -r -o 2019-12803-6/sketch -I 2019-12803-6 -C 2019-12803-6/R1.fq.gz 2019-12803-6/R1.fq.gz