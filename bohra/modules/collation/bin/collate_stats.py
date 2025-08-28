#!/usr/bin/env python3
import sys, json, datetime, subprocess, pathlib,pandas


def get_vals_seqtk(f):

    with open(f, 'r') as inputfile:
        lines = inputfile.read().split('\n')
        qscore = lines[1].split('\t')[3]
        bases = lines[1].split('\t')[0]
        gc= float(lines[1].split('\t')[1]) + float(lines[1].split('\t')[2])
        try:
            return round(float(qscore),1)
        except ValueError:
            print(f"Something has gone horribly wrong {qscore} can not be parsed to a float!!")
            return 0.0
        
def sum_stats(f,cmd,col):

    p = subprocess.run(f"cat {f} | datamash {cmd} --header-in {col}", shell = True, capture_output = True, encoding= 'utf-8')

    return p.stdout

def get_vals_kit(f):

    reads = sum_stats(f, 'sum', 4)
    minval = sum_stats(f, 'mean', 6)
    avgval = sum_stats(f, 'mean', 7)
    maxval = sum_stats(f, 'mean', 8)

    return int(reads.strip()), float(minval.strip()), float(avgval.strip()), float(maxval.strip())

def get_length(genome_size):

    with open(genome_size, 'r') as f:
        x = f.read().strip().split('\n')[0].strip()
        print(x)
        try:
            return int(x)
        except:
            print(f"Something has gone horribly wrong {x} can not be parsed to an int!!")
            return 0
    
def get_dpth(genome_size,bases):
    length = get_length(genome_size)
    try:
        dpth = int(bases)/length
        return round(dpth, 1),int(length)
    except:
        print(f"Something has gone horribly wrong with the genome size {genome_size} {bases} can not be parsed to an int!!")
        return 0.0, int(length)
    
def get_gcs(path):

    try:
        gcs = pandas.read_csv(path, sep = '\t')
        return gcs[gcs.columns[0]].values[0]
    except Exception as e:
        print(f"Something has gone wrong with the GC content file {path} {e}")
        return 0.0

tab = pandas.read_csv(sys.argv[2], sep = '\t')
tab['Isolate'] = sys.argv[1]
print(tab)
gcs = get_gcs(sys.argv[3])
print(gcs)
print(sys.argv[6])
tab = tab.rename(columns = {'num_seqs': 'Reads', 'sum_len': 'Yield','min_len':'Min len','max_len':'Max len', 'avg_len':'Avg len', 'Q30(%)': 'Qscore'})
dpth,size = get_dpth(genome_size = sys.argv[4], bases = tab['Yield'].values[0])
tab['Depth'] = dpth
tab["Genome size"] = size
tab['GC'] = gcs
tab['is_control'] = True if "control" in sys.argv[7] else False
tab['filesize'] = "<20000000" if sys.argv[8] == "FAIL_READ_FILE_TOO_SMALL" else ">20000000"
tab['Qscore'] = get_vals_seqtk(sys.argv[5])

tab = tab[['Isolate','Reads','Yield','GC','Min len','Avg len','Max len','Qscore',"Genome size", 'Depth',"is_control","filesize"]]
tab.to_csv('read_assessment.txt', sep = '\t', index = False)





    


