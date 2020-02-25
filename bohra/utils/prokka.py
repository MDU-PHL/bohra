import toml, pathlib, subprocess, sys

def generate_prokka_cmd(isolate, assembly):
    
    cm f"prokka --outdir {isolate} --prefix {isolate} --mincontiglen 500 --notrna --fast --force {assembly} --cpus 1"

    return cmd

def generate_rm_cmd(isolate):

    cmd = f"rm {isolate}/*.err {isolate}/*.faa {isolate}/*.ffn {isolate}/*.fsa {isolate}/*.sqn {isolate}/*.tbl {isolate}/*.tsv"
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
    
def main(inputs, isolate):
    
    # set up data dict
    data = {}
    data[isolate] = {}
    data[isolate]['prokka'] = {}
    # run kraken
    cmd = generate_prokka_cmd(isolate = isolate, assembly = inputs)
    p = run_cmd(cmd)
    if p == 0:
        rm_cmd = generate_rm_cmd(isolate = isolate)
        r = generate_rm_cmd(run_cmd)
        if r == 0:
            data[isolate]['prokka'] = {}
            data[isolate]['prokka']['gff'] = f'{sample}/{sample}.gff'
            data[isolate]['prokka']['txt'] = f'{sample}/{sample}.txt'
            [f'{sample}.gff', f'{sample}.txt']
            write_toml(data = data, output= f'{isolate}/prokka.toml')



if __name__ == '__main__':
    
    main(inputs = f"{sys.argv[1]}", isolate = f"{sys.argv[2]}")
    

