import toml, pathlib, subprocess, sys, pandas
from snakemake import shell
from Bio import SeqIO
                
                
def check_snippy(data, isolate, minaln, aln):
    
    p = pathlib.Path(aln)
    fasta = p.open()
    length = 0
    nocov = 0
    lowcov = 0
    het = 0
    unaln = 0
    for i in SeqIO.parse(fasta,'fasta'): # use BioPython to determine percent alignment
        length = length + len(i.seq)
        nocov = nocov + i.seq.count('-')
        lowcov = lowcov + i.seq.count('N')
        het = het + i.seq.count('n')
    unaln = unaln + nocov + lowcov + het
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
        print(f'Isolate {isolate} will be included in core analysis')
    else:
        data[isolate]['qc_snippy']['Quality'] = f'Alignment < {float(minaln)} - isolate will not be included in core'
        print(f'Isolate {isolate} did not meet the {minaln} threshold for alignment and will not be included in the core analysis.')
    return data
                                        

def open_toml(tml):

    data = toml.load(tml)
    return data

def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)
    
def main(inputs, isolate, output, minaln):
    

    s = open_toml(tml = inputs)
    print(f'Checking if snippy was performed for isolate {isolate}.')
    rs = s[isolate]['snippy']['done']
    data = {}
    data[isolate] = {}
    data[isolate]['qc_snippy'] = {}
    data[isolate]['qc_snippy']['run_snippy'] = rs
    
    if rs == 'Yes':
        print(f'snp calling was performed on isolate {isolate}, further checks will now be performed.')
        data = check_snippy(minaln = minaln,data = data, aln = s[isolate]['snippy']['alignment'], isolate = isolate)
    else:
        data[isolate]['qc_snippy']['Quality'] = 'FAILED sequence QC will not be included in further analysis.'
    print(f'Alignment assessment has been performed for isolate {isolate} and is now being written to toml.')
    write_toml(data = data, output = f"{isolate}/snippy_qc.toml") 

   
# {input} {wildcards.sample} {output} {params.minaln}
inputs = snakemake.input
isolate = snakemake.wildcards.sample
output = snakemake.output
minaln=snakemake.params.minaln

main(inputs = inputs, isolate = isolate, output = output, minaln = minaln)
