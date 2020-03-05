import toml, pathlib, subprocess, sys, pandas

def generate_abritamr_cmd(input_file, isolate, wd):
    
    w = pathlib.Path(wd)
    cmd = f"abriTAMR -c {input_file} -pfx {isolate} -w {w}"
    return cmd

def generate_dict(isolate):

    df = pandas.read_csv(f"{isolate}/summary_matches.csv")
    return df.to_dict(orient = 'records')

def run_cmd(cmd):
    
    p = subprocess.run(cmd, shell = True, capture_output=True, encoding = 'utf-8')
    return p.stderr

def open_toml(tml):

    data = toml.load(tml)

    return data

def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)
    
def main(inputs, isolate, seqdata, wd):
    
    seqdata = open_toml(seqdata)
    data = {}
    data[isolate] = {}
    data[isolate]['resistome'] = {}
    if seqdata[isolate]['seqdata']['data']['Quality'] == 'PASS':
    # set up data dict
        print('doing something')
        tml = open_toml(inputs)
        contigs = f"{isolate}/contigs.fa"
        cmd = generate_abritamr_cmd(input_file = contigs, isolate = isolate, wd = wd)
        print(cmd)
        p = run_cmd(cmd)
        data[isolate]['resistome']['tool'] = 'abritamr'
        
        if f"pipeline successfully completed" in p:
            data[isolate]['resistome']['data'] = generate_dict(isolate)
            data[isolate]['resistome']['done'] = True
        else:
            data[isolate]['resistome']['done'] = False
    else:
        data[isolate]['resistome']['done'] = False
        data[isolate]['resistome']['tool'] = f"AMR not performed - failed QC"
    
    write_toml(data = data, output= f'{isolate}/resistome.toml')



if __name__ == '__main__':
    
    main(inputs = f"{sys.argv[1]}", isolate = f"{sys.argv[2]}", seqdata = f"{sys.argv[3]}", wd = f"{sys.argv[4]}")
    

