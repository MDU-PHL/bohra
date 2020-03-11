import toml, pathlib, subprocess, sys, pandas
from snakemake import shell
from Bio import SeqIO
                
                
def check_snippy(data, isolate, minaln, aln):
    
    p = pathlib.Path(aln)
    fasta = p.open()
    for i in SeqIO.parse(fasta,'fasta'): # use BioPython to determine percent alignment
        length = len(i.seq)
        nocov = i.seq.count('-')
        lowcov = i.seq.count('N')
        het = i.seq.count('n')
        unaln = nocov + lowcov + het
        perc_aln = 100*(length - unaln) / length
        data[isolate]['qc_snippy']['length'] = length
        data[isolate]['qc_snippy']['nocov'] = nocov
        data[isolate]['qc_snippy']['lowcov'] = lowcov
        data[isolate]['qc_snippy']['het'] = het
        data[isolate]['qc_snippy']['unaln'] = unaln
        data[isolate]['qc_snippy']['pecr_aln'] = perc_aln
        # if the percent alignement is greater than the min alignment
        if perc_aln > float(minaln):
            data[isolate]['qc_snippy']['Quality'] = 'PASS'
        else:
            data[isolate]['qc_snippy']['Quality'] = f'Alignment < {float(minaln)} - isolate will not be included in core'
        return data
                                        

def open_toml(tml):

    data = toml.load(tml)
    return data

def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)
    
def main(inputs, isolate, output, minaln):
    

    s = open_toml(tml = inputs)
    
    rs = s[isolate]['snippy']['done']
    data = {}
    data[isolate] = {}
    data[isolate]['qc_snippy'] = {}
    data[isolate]['qc_snippy']['run_snippy'] = rs
    
    if rs == 'Yes':
        data = check_snippy(minaln = minaln,data = data, aln = s[isolate]['snippy']['alignment'], isolate = isolate)
    else:
        data[isolate]['qc_snippy']['Quality'] = 'FAILED sequence QC will not be included in further analysis.'
    
    write_toml(data = data, output = f"{isolate}/snippy_qc.toml") 

   
# {input} {wildcards.sample} {output} {params.minaln}
inputs = snakemake.input
isolate = snakemake.wildcards.sample
output = snakemake.output
minaln=snakemake.params.minaln

main(inputs = inputs, isolate = isolate, output = output, minaln = minaln)
