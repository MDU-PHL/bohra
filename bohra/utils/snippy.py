import toml, pathlib, subprocess, sys


def generate_snippy_cmd(r1, r2, isolate, reference, threads):
    
    cmd = f"snippy --outdir {isolate} --ref {reference} --R1 {r1} --R2 {r2} --force --cpus {threads}"
    print(cmd)
    return cmd

def run_cmd(cmd):
    
    p = subprocess.run(cmd, shell = True, capture_output=True, encoding = 'utf-8')
    return p.returncode

def get_reads(inputs, isolate):

    s = open_toml(inputs)
    r1 = s[isolate]['seqdata']['R1']
    r2 = s[isolate]['seqdata']['R2']
    return r1,r2

def get_quality(inputs, isolate):

    s = open_toml(inputs)
    print(s)
    if s[isolate]['seqdata']['data']['Quality'] == 'PASS':
        return 'Yes'
    else:
        return 'No'


def open_toml(tml):

    data = toml.load(tml)

    return data

def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)
    
def main(inputs, isolate, output, reference, threads):
    
    print(inputs)
    # print(threads)
    r1,r2 = get_reads(inputs = inputs, isolate = isolate)
    run_snippy = get_quality(inputs = inputs, isolate = isolate)
    print(run_snippy)
    data = {}
    data[isolate] = {}
    data[isolate]['snippy'] = {}
    data[isolate]['snippy']['reference'] = reference
    data[isolate]['snippy']['R1'] = r1
    data[isolate]['snippy']['R2'] = r2
    data[isolate]['snippy']['run_snippy'] = run_snippy

    if run_snippy == 'Yes':
        cmd = generate_snippy_cmd(r1 = r1, r2=r2, isolate = isolate, reference = reference, threads = threads)
        # print(cmd)
        p = run_cmd(cmd)
        # print(p)
        if p == 0:
            data[isolate]['snippy']['alignment'] = f"{isolate}/snps.aligned.fa"
            data[isolate]['snippy']['vcf'] = f"{isolate}/snps.vcf"
            data[isolate]['snippy']['cmd'] = cmd
    write_toml(data = data, output = f"{isolate}/snippy.toml") 
    

if __name__ == '__main__':
    
    main(inputs = f"{sys.argv[1]}", isolate = f"{sys.argv[2]}", output = f"{sys.argv[3]}",reference = f"{sys.argv[4]}", threads =f"{sys.argv[5]}" )
    



# mash triangle -C *.msh

# mash sketch -m 5 -s 10000 -r -o 2019-12803-6/sketch -I 2019-12803-6 -C 2019-12803-6/R1.fq.gz 2019-12803-6/R1.fq.gz