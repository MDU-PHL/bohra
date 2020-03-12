import toml, pathlib, subprocess, sys, pandas, datetime
from snakemake import shell

def generate_iqtree_cmd(aln, script_path):

    cmd = f"bash {script_path}/iqtree_generator.sh ref.fa {aln} core 20"								
    return cmd

def generate_delete_cmd():
    cmd = f"rm *.ckp.gz *.contree *.bionj"
    return cmd

def run_cmd(cmd):

    print(f"Running : {cmd}.")
    p = subprocess.run(cmd, shell = True, capture_output=True, encoding = 'utf-8')
    return p.stdout

def open_toml(tml):

    data = toml.load(tml)

    return data

def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)
    
def main(inputs, ref, idx, script_path):
    
    t = open_toml(inputs)
    if t['gubbins']['run_gubbins'] == 'Yes':
        aln = 'gubbins.aln'
    else:
        aln = 'core.aln'
    print(f"Generating iqtree command.")
    cmd = generate_iqtree_cmd(aln = aln, script_path = script_path)
    i = run_cmd(cmd)
    i = f"{i.strip()} -redo"
    iqtree = run_cmd(i)
    # print(iqtree)
    dl = generate_delete_cmd()
    rm = run_cmd(dl)
    data = {}
    data['iqtree'] = {}
    data['iqtree']['command'] = cmd
    data['iqtree']['input_file'] = aln
    data['iqtree']['file'] = 'core.treefile'

    write_toml(data = data, output = "iqtree.toml")

inputs = snakemake.input.gubbins
ref = snakemake.params.ref
idx = snakemake.params.idx
script_path = snakemake.params.script_path

main(inputs = inputs, ref = ref, idx = idx, script_path = script_path)

