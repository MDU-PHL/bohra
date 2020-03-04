import toml, pathlib, subprocess, sys, pandas

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
    print(data)
    return data

def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)
    
def main(inputs, isolate, output, minaln):
    

    s = open_toml(tml = inputs)
    print(s)
    rs = s[isolate]['snippy']['run_snippy']
    data = {}
    data[isolate] = {}
    data[isolate]['qc_snippy'] = {}
    data[isolate]['qc_snippy']['run_snippy'] = rs
    if rs:
        data = check_snippy(minaln = minaln,data = data, aln = s[isolate]['snippy']['alignment'], isolate = isolate)
    else:
        data[isolate]['qc_snippy']['Quality'] = 'FAILED sequence QC will not be included in further analysis.'
    
    write_toml(data = data, output = f"{isolate}/snippy_qc.toml") 

   

if __name__ == '__main__':
    
    main(inputs = f"{sys.argv[1]}", isolate = f"{sys.argv[2]}", output = f"{sys.argv[3]}",minaln = f"{sys.argv[4]}")
    



# mash triangle -C *.msh

# mash sketch -m 5 -s 10000 -r -o 2019-12803-6/sketch -I 2019-12803-6 -C 2019-12803-6/R1.fq.gz 2019-12803-6/R1.fq.gz