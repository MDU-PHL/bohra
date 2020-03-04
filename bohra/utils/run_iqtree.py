import toml, pathlib, subprocess, sys, pandas, datetime

def generate_iqtree_cmd(script_path, aln):

    cmd = f"bash {script_path}/iqtree_generator.sh ref.fa {aln} core 20"								
    return cmd

def generate_delete_cmd():
    cmd = f"rm *.ckp.gz *.contree *.bionj"
    return cmd

def run_cmd(cmd):

    p = subprocess.run(cmd, shell = True, capture_output=True, encoding = 'utf-8')
    print(p)
    return p.stdout

def open_toml(tml):

    data = toml.load(tml)

    return data

def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)
    
def main(inputs, ref, idx, script_path):
    
    t = open_toml(inputs)
    print(t)
    if t['gubbins']['run_gubbins']:
        aln = 'gubbins.aln'
    else:
        aln = 'core.aln'
    cmd = generate_iqtree_cmd(aln = aln, script_path = script_path)
    i = run_cmd(cmd)
    i = f"{i.strip()} -redo"
    print(i)
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

if __name__ == '__main__':
    
    main(inputs = f"{sys.argv[1]}", ref = f"{sys.argv[2]}", idx = f"{sys.argv[3]}",script_path = f"{sys.argv[4]}")
    



# mash triangle -C *.msh

# mash sketch -m 5 -s 10000 -r -o 2019-12803-6/sketch -I 2019-12803-6 -C 2019-12803-6/R1.fq.gz 2019-12803-6/R1.fq.gz