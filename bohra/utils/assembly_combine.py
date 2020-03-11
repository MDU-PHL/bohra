import toml, pathlib, subprocess, sys, pandas, numpy

from snakemake import shell
import pandas, pathlib


def assembly_combine(assemblies):

    ac = pandas.DataFrame()
    
# colnames = ['Name','bp','# Contigs','Ns','# Gaps','Min Contig size','Max Contig size','Avg Contig size','N50']
# {'Name': isolate, 'bp':length, '# Contigs':no, 'Ns':Ns, '# Gaps':gaps, 'Min Contig size':min_len, 'Max Contig size':max_len,  'Avg Contig size':avg_len, 'N50':N50}
    for a in assemblies:
        d = open_toml(a)
        isolate=list(d.keys())[0]
        df = pandas.DataFrame(d[isolate]['assembly_stats'], index = [0])
        if ac.empty:
            ac = df
        else:
            ac = df.append(ac, sort = True)
    # print(ac)
    return ac


def combine(prokka, assembly_stats):
# 
    
    gff = pandas.DataFrame()
    
    for p in prokka:
        tml = open_toml(p)
        isolate=list(tml.keys())[0]
        if tml[isolate]['prokka']['done'] and 'txt' in tml[isolate]['prokka']:
            g = pathlib.Path(tml[isolate]['prokka']['txt'])
            df = pandas.read_csv(g, sep = ':', header = None, names = ['cond', f"{g.parts[0]}"])
            # print(df)
            if gff.empty:
                    gff = df
            else:
                    gff = gff.merge(df, how = 'outer')
    gff = gff[gff['cond'].isin(['CDS', 'rRNA'])]
    gff = gff.T
    gff.columns = gff.iloc[0]
    gff = gff.iloc[1:]  
    assembly = assembly_stats.merge(gff, left_on = ['Name'], right_on= gff.index, how = 'outer')

    m = list(assembly[assembly['# Contigs'] != '-']['# Contigs'])  
    mn = numpy.mean(m)
    s =(2* numpy.std(m))
    cut = mn + s
    assembly['Quality'] = numpy.where(assembly['# Contigs'].replace('-',0) <= cut, assembly['Quality'], f'WARNING - assembly outlier (> {round(cut,2)} contigs)')
    assembly = assembly.rename(columns={'Name':'Isolate'})
    
    assembly = assembly[['Isolate', 'bp','# Contigs','Ns','# Gaps','Min Contig size','Max Contig size','Avg Contig size','N50', 'CDS','rRNA', 'Quality' ]]

    assembly.to_csv(f"assembly.tab", sep = '\t', index = False)

    return assembly.to_dict(orient = 'records')

def open_toml(tml):

    data = toml.load(tml)

    return data

def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)
    
def main(prokka):
    
    assemblies = assembly_combine(assemblies = prokka)
    assembly_stat_dict = combine(prokka = prokka, assembly_stats = assemblies) 
    # set up data dict
    data = {}
    data['assembly_stat'] = assembly_stat_dict
    write_toml(data = data, output= f'assembly.toml')


prokka = snakemake.input.prokka

main(prokka = prokka)
    

