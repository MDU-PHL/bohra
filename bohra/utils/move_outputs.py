import toml, pathlib, subprocess, sys
from snakemake import shell

def generate_mv_cmd(iostring):

    cmd = f"mv {iostring} && mv *.tab report/"
    return cmd

def get_io(inputs):

    i = f"{inputs}"
    o = f"{pathlib.Path('report', inputs)}" if f"{inputs}" != 'report.html' else f"{pathlib.Path('report', 'index.html')}"
    return f"{i} {o}"

def run_cmd(cmd):
    
    p = subprocess.run(cmd, shell = True, capture_output=True, encoding = 'utf-8')
    return p.returncode

    
def main(inputs):
    
    for i in inputs:
        iostring = get_io(i)
        cmd = generate_mv_cmd(iostring)
        print(f"Moving out files : {cmd}")
        p = run_cmd(cmd)
    


inputs = snakemake.input    
main(inputs = inputs)
    

