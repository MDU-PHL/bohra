from Bio import SeqIO
import pathlib, toml, pandas, sys
from snakemake import shell


def open_toml(tml):

    data = toml.load(tml)

    return data

def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)


def fa(inputs, isolate, min_size =500):
    p = pathlib.Path(f"{inputs}")
    min_len = 1.0e+10
    max_len = 0
    contig_length = []
    Ns = 0
    gaps = 0
    contigs_less_than_min = 0
    min_size = int(min_size)
    no = 0
    for i in SeqIO.parse(f"{p}", 'fasta'):
            contig_length.append(len(i.seq))
            if len(i.seq) < min_size:
                    contigs_less_than_min += 1
            else:   
                    min_len = len(i.seq) if len(i.seq) < min_len else min_len
                    max_len = len(i.seq) if len(i.seq) > max_len else max_len
                    Ns = Ns + i.seq.lower().count('n')
                    gaps = gaps + i.seq.count('-')
                    no += 1
    length = sum([c for c in contig_length if c > min_size])
    len_sorted = sorted([c for c in contig_length if c > min_size])
    avg_len = int(length/ len([c for c in contig_length if c > min_size]))
    cum = 0 
    N50 = 0
    for c in len_sorted:
            # print(c)
            cum += c
            # print(cum)
            if cum >= length/2:
                    N50 = c
                    break

    data = {'Name': isolate, 'bp':length, '# Contigs':no, 'Ns':Ns, '# Gaps':gaps, 'Min Contig size':min_len, 'Max Contig size':max_len,  'Avg Contig size':avg_len, 'N50':N50, 'Quality': 'PASS'}
    return(data)

def main(inputs,isolate):

    assembly_data = open_toml(tml = inputs)

    if assembly_data[isolate]['assembly']['done'] == 'Yes':
        # print('in true')
        print(f"Calculating assembly statistics.")
        assembly = f"{pathlib.Path(isolate, 'contigs.fa')}"
        d = fa(inputs=assembly,isolate= isolate)
        data = {}
        data[isolate]={}
        data[isolate]['assembly_stats'] = d
        
    else:
        data = {}
        data[isolate]={}
        data[isolate]['assembly_stats'] = {'Name': isolate, 'bp':'-', '# Contigs':'-', 'Ns':'-', '# Gaps':'-', 'Min Contig size':'-', 'Max Contig size':'-',  'Avg Contig size':'-', 'N50':'-', 'Quality': 'ASSEMBLY NOT PERFORMED please see Sequence data tab'}
    write_toml(data= data, output = f'{isolate}/assembly_stats.toml')
        
 
inputs = snakemake.input
isolate = snakemake.wildcards.sample

main(inputs = inputs, isolate = isolate)
    


