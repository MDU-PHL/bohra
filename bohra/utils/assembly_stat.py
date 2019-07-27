from Bio import SeqIO
import pathlib
import pandas, argparse



def fa(path, min_size, full_path):
    p = pathlib.Path(f"{path}")
    isolate = f"{p}" if full_path else p.name.split('.')[0]
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

def get_fa_stat(assemblies, minsize, full_path):
    colnames = ['Name','bp','# Contigs','Ns','# Gaps','Min Contig size','Max Contig size','Avg Contig size','N50']
    # print(colnames)
    print('\t'.join(colnames))
    for a in assemblies:
        data = fa(path=a, min_size = minsize, full_path = full_path)
        # print(data)
        print('\t'.join([f"{data[x]}" for x in colnames]))
        


def set_parsers():
    parser = argparse.ArgumentParser(description='Python version of torstyvers fa',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('assemblies', help = 'white space separated list of assemblies', nargs='+')
    parser.add_argument('-m', '--minsize', help = 'min contig size to be included in calculations', default = 500)
    parser.add_argument('-f', '--full_path',action = 'store_true', help = 'Set if you would like to use the full path as the isolate identifer - default to the name of the fasta file')
    


    args = parser.parse_args()
    return(args)

def main():

    args = set_parsers()
    if vars(args) == {}:
        parser.print_help(sys.stderr)
    else:
        assemblies = args.assemblies
        full_path = args.full_path
        minsize = args.minsize
        get_fa_stat(assemblies, minsize, full_path)


if __name__ == '__main__':
    main()


