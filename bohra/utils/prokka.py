import toml, pathlib, subprocess, sys

def generate_prokka_cmd(isolate, assembly):
    
    cmd =  f"prokka --outdir {isolate} --prefix {isolate} --mincontiglen 500 --notrna --fast --force {assembly} --cpus 4"

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
    
def main(inputs, isolate, seqdata):
    
    print(inputs)
    # set up data dict
    data = open_toml(inputs)
    # data[isolate] = {}
    data[isolate]['prokka'] = {}
    # run kraken
    assembly = f"{isolate}/contigs.fa"
    seqdata = open_toml(seqdata)
    cmd = generate_prokka_cmd(isolate = isolate, assembly = assembly)
    if seqdata[isolate]['seqdata']['data']['Quality'] == 'PASS':
        print(cmd)
        p = run_cmd(cmd)
        if p == 0:
            rm_cmd = generate_rm_cmd(isolate = isolate)
            print(rm_cmd)
            r = generate_rm_cmd(run_cmd)
            data[isolate]['prokka']['done'] = True
            data[isolate]['prokka']['gff'] = f'{isolate}/{isolate}.gff'
            data[isolate]['prokka']['txt'] = f'{isolate}/{isolate}.txt'
    else:
        data[isolate]['prokka']['done'] = False

    
    write_toml(data = data, output= f'{isolate}/prokka.toml')



if __name__ == '__main__':
    
    main(inputs = f"{sys.argv[1]}", isolate = f"{sys.argv[2]}", seqdata = f"{sys.argv[3]}")
    

