from Bio import SeqIO
import pathlib, toml, pandas



def open_toml(tml):

    data = toml.load(tml)

    return data

def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)


def fa(inputs, isolate, min_size =500):
    p = pathlib.Path(f"{path}")
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

    data = {'Name': isolate, 'bp':length, '# Contigs':no, 'Ns':Ns, '# Gaps':gaps, 'Min Contig size':min_len, 'Max Contig size':max_len,  'Avg Contig size':avg_len, 'N50':N50}
    return(data)

def get_fa_stat(assembly,isolate):
#     colnames = ['Name','bp','# Contigs','Ns','# Gaps','Min Contig size','Max Contig size','Avg Contig size','N50']
    # print(colnames)
    # print('\t'.join(colnames))
    d = fa(path=assembly, isolate)
    data = {}
    data[isolate]={}
    data[isolate]['assembly_stats'] = d
    write_toml(data= data, output = f'{isolate}/assembly_stats.toml')
        # print(data)
        # print('\t'.join([f"{data[x]}" for x in colnames]))
        
def main():

    main(inputs = f"{sys.argv[1]}", minsize = f"{sys.argv[2]}", isolate = f"{sys.argv[3]}",  prefill =  f"{sys.argv[4]}")
    


if __name__ == '__main__':
    main()


