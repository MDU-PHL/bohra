import toml, pathlib, subprocess, sys, pandas
from snakemake import shell

def open_toml(tml):
    
    data = toml.load(tml)
    
    return data
    
def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)
    
def main(inputs):
    data = {}
    tab = pandas.DataFrame()
    for i in inputs:
        # print(i)
        d = open_toml(i)
        # print(list(d.keys()))
        isolate = list(d.keys())[0]
        df = pandas.DataFrame(data = d[isolate]['seqdata']['data'], index = [0])
        df['Isolate'] = isolate
        # print(df)
        if tab.empty:
            tab = df
        else:
            tab = tab.append(df)
    tab = tab.rename(columns = {'min_len': 'Min len', 'max_len': 'Max len', 'avg_len': 'Avg len', 'avgQ': 'Avg Qual', 'Estimated_coverage': 'Estimated depth', 'bases':'Yield'})
    tab = tab[['Isolate','Reads','Yield','GC content','Min len','Avg len','Max len','Avg Qual','Estimated depth', 'Quality']]
    tab.to_csv(f"seqdata.tab", sep = '\t', index = False)
    data['seqdata'] = tab.to_dict(orient = 'records')

    write_toml(data = data, output= "seqdata.toml")
    

inputs = snakemake.input

main(inputs = inputs)
# if __name__ == '__main__':
    
#     main(inputs = sys.argv[1:])
    




# mash triangle -C *.msh

# mash sketch -m 5 -s 10000 -r -o 2019-12803-6/sketch -I 2019-12803-6 -C 2019-12803-6/R1.fq.gz 2019-12803-6/R1.fq.gz