import toml, pathlib, subprocess, sys, pandas, json


def get_gffs(inputs):

    l = []
    for i in inputs:
        tml = open_toml(i)
        isolate = tml.keys()[0]
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
    
    roary = run_cmd(generate_roary_cmd(gffs))
    if roary == 0:
        generate_mv_cmd()
        generate_rm_cmd()

        data['roary']['done'] = True
        data['roary']['csv'] = "roary/gene_presence_absence.csv"
        data['roary']['txt'] "roary/summary_statistics.txt"
    
    
    write_toml(data = data, output= f'roary.toml')



if __name__ == '__main__':
    
    main(inputs = sys.argv[1])
    

