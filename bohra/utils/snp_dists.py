import toml, pathlib, subprocess, sys, pandas, datetime

def generate_dists_cmd(aln):
    
    cmd = f"snp-dists {aln} > distances.tab"								
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
    
def main(inputs):
    
    t = open_toml(inputs)
    
    if t['gubbins']['run_gubbins']:
        aln = 'gubbins.aln'
    else:
        aln = 'core.aln'
    
    dists= run_cmd(generate_dists_cmd(aln = aln))
    if dists == 0:
        data = {}
        data['snp_dists'] = {}
        data['snp_dists']['input_file'] = aln
        data['snp_dists']['file'] = 'distances.tab'


        
        write_toml(data = data, output = "distances.toml")

if __name__ == '__main__':
    
    main(inputs = f"{sys.argv[1]}")
    



# mash triangle -C *.msh

# mash sketch -m 5 -s 10000 -r -o 2019-12803-6/sketch -I 2019-12803-6 -C 2019-12803-6/R1.fq.gz 2019-12803-6/R1.fq.gz