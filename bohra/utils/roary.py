import toml, pathlib, subprocess, sys, pandas, json
from snakemake import shell

def get_gffs(inputs):

    l = []
    for i in inputs:
        tml = open_toml(i)
        isolate = list(tml.keys())[0]
        if tml[isolate]['prokka']['done']:
            gff = tml[isolate]['prokka']['gff']
            l.append(gff)
    return ' '.join(l)

def generate_roary_cmd(gffs):

    cmd = f"nice roary -p 36 -f roary {gffs}"
    return cmd

def generate_mv_cmd():
    
    cmd = f"mv roary_*/* roary"
    return cmd

def generate_rm_cmd():
    
    cmd = f"rm -r roary_*"
    return cmd

def run_cmd(cmd):
    
    print(f"Running : {cmd}.")
    p = subprocess.run(cmd, shell = True, capture_output=True, encoding = 'utf-8')
    return p.returncode

def open_toml(tml):

    data = toml.load(tml)

    return data

def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)
    
def main(inputs):
    
    gffs = get_gffs(inputs)
    data = {}
    data['roary'] = {}
    print(f"Roary will now be performed.")
    roary = run_cmd(generate_roary_cmd(gffs))
    if roary == 0:
        print(f"Roary successfully completed.")
        run_cmd(generate_mv_cmd())
        run_cmd(generate_rm_cmd())

        data['roary']['done'] = True
        data['roary']['csv'] = "roary/gene_presence_absence.csv"
        data['roary']['txt']="roary/summary_statistics.txt"
    else:
        print('Something went wrong - roary did not complete.')
        data['roary']['done'] = False
    
    write_toml(data = data, output= f'roary.toml')

inputs = snakemake.input.prokka

main(inputs = inputs)
    

