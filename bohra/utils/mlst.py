import toml, pathlib, subprocess, sys, pandas, json

def generate_mlst_cmd(assembly, isolate):

    cmd = f"mlst --label {isolate} --json {isolate}.json {assembly}"
    return cmd

def extract_mlst(isolate):

    jsn = json.load(open(f"{isolate}/{isolate}.json"))
    isolate = jsn[0]['id']
    alleles = []
    ax = zip(jsn[0]['alleles'].keys(), jsn[0]['alleles'].values())
    for a in ax:
        alleles.append(f"{a[0]}({a[1]})")
    alleles = sorted(alleles)

    d = {
        'Isolate':isolate,
        'alleles': alleles,
        'len_alleles': len(alleles),
        'scheme':jsn[0]['scheme'],
        'ST':jsn[0]['sequence_type']
    }
    return d

def generate_rm_cmd(isolate):

    cmd = f"rm {isolate}/{isolate.json}"
    return df.to_dict(orient = 'records')

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
    
    data = {}
    seqdata = open_toml(seqdata)
    if seqdata[isolate]['seqdata']['data']['Quality'] == 'PASS'
    # set up data dict
        cmd = generate_abritamr_cmd(input_file = inputs, isolate = isolate)
        p = run_cmd(cmd)
        if p == 0:
            data[isolate]['mlst']['done'] = True
            data[isolate]['mlst']['data'] = extract_mlst(isolate)
    else:
            data[isolate]['mlst']['done'] = False
            data[isolate]['mlst']['data'] = f"MLST not performed - failed QC"
            
    write_toml(data = data, output= f'{isolate}/mlst.toml')



if __name__ == '__main__':
    
    main(inputs = f"{sys.argv[1]}", isolate = f"{sys.argv[2]}", seqdata = f"{sys.argv[3]}")
    

