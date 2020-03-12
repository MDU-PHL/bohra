import toml, pathlib, subprocess, sys, pandas
from snakemake import shell

def generate_triangle_cmd():
    
    cmd = f"mash triangle -C */sketch.msh ref.fa > preview_distances.tab"
    return cmd

def generate_tree_cmd():

    cmd = f"quicktree -in m -out t preview_distances.tab  | nw_order -c n - > preview.newick"
    return cmd

def run_cmd(cmd):
    
    print(f"Running : {cmd}")
    p = subprocess.run(cmd, shell = True, capture_output=True, encoding = 'utf-8')
    print(p)
    return p.returncode

def get_isolates(tab):
    data = open(tab).readlines()
    isolates = ['mash-dist']
    for d in range(1,len(data)):
        l = data[d].split()
        isolates.append(l[0])
    return isolates

def clean_dist_mat(tab):
    isolates = get_isolates(tab)
    df = pandas.read_csv(tab, names = isolates, sep = '\t', header = None, skiprows = [0])
    df = df.fillna('-')
    df.to_csv(tab, sep = '\t', index = False)

def open_toml(tml):
    # print(tml)
    data = toml.load(tml)

    return data
    
def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)
    
def main(inputs):
    
    triangle = generate_triangle_cmd()
    # print(triangle)
    p = run_cmd(triangle)
    # print(p)
    print(f"Generating preview files.")
    if p == 0:
        tree = generate_tree_cmd()
        # print(tree)
        q = run_cmd(tree)
        if q == 0:
            data = {}
            data['preview'] = {}
            data['preview']['tree'] = 'preview.newick'
            data['preview']['distances'] = 'preview_distances.tab'
            clean_dist_mat(data['preview']['distances'])
            write_toml(data = data, output = 'preview.toml')
    
# reference = snakemake.params.reference
inputs = snakemake.input

main(inputs = inputs)

# if __name__ == '__main__':
    
#     main(inputs = sys.argv[1:])
    



# mash triangle -C *.msh

# mash sketch -m 5 -s 10000 -r -o 2019-12803-6/sketch -I 2019-12803-6 -C 2019-12803-6/R1.fq.gz 2019-12803-6/R1.fq.gz