import toml, pathlib, subprocess, sys, pandas, datetime

def generate_snippy_clean_cmd():

    cmd = f"snippy-clean_full_aln core.full.aln  > clean.full.aln"								
    return cmd

def generate_gubbins_cmd():

    cmd = f"run_gubbins.py -c 36  --prefix core clean.full.aln"								
    return cmd

def generate_snp_sites_cmd():
    
    cmd = f"run_gubbins.py -c 36  --prefix core gubbins.aln"								
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
    
def main(inputs, gubbins):
    data = {}
    data['gubbins'] = {}
    # print(gubbins)
    if gubbins != 'False':
        print('inside true')
        clean = run_cmd(generate_snippy_clean_cmd())
        if clean == 0:
            gub = run_cmd(generate_gubbins_cmd())
            if gub == 0:
                sites = run_cmd(generate_snp_sites_cmd())
                if sites == 0:
                    data['gubbins']['run_gubbins'] = True
                    data['gubbins']['aln'] = 'gubbins.aln'
    else:
        data['gubbins']['run_gubbins'] = False

        
    write_toml(data = data, output = "gubbins.toml")

if __name__ == '__main__':
    
    main(inputs = f"{sys.argv[1]}", gubbins = sys.argv[2])
    



# mash triangle -C *.msh

# mash sketch -m 5 -s 10000 -r -o 2019-12803-6/sketch -I 2019-12803-6 -C 2019-12803-6/R1.fq.gz 2019-12803-6/R1.fq.gz