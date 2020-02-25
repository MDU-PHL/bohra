import toml, pathlib, subprocess, sys, pandas


import pandas, pathlib


def assembly_combine(assemblies):

    ac = pandas.DataFrame()
    
# colnames = ['Name','bp','# Contigs','Ns','# Gaps','Min Contig size','Max Contig size','Avg Contig size','N50']
# {'Name': isolate, 'bp':length, '# Contigs':no, 'Ns':Ns, '# Gaps':gaps, 'Min Contig size':min_len, 'Max Contig size':max_len,  'Avg Contig size':avg_len, 'N50':N50}
    for a in assemblies:
        d = open_toml(a)
        isolate=d.keys()[0]
        df = pandas.DataFrame(d =d[isolate], index = [0])
        df['Quality'] = 'PASS'
        if not tml[isolate]['resistome']['done']:
            df['Quality '] = 'NOT INCLUDED - FAIL Sequence QC'
        if ac.empty:
            ac = df
        else:
            ac = df.append(ac, sort = True)
    
    return ac


def combine(prokka, assembly_stats)
# 
    
    gff = pandas.DataFrame()
    
    for p in prokka:
        tml = open_toml(p)
        isolate=tml.keys()[0]
        g = pathlib.Path(tml[isolate]['prokka']['gff'])
        df = pandas.read_csv(g, sep = ':', header = None, names = ['cond', f"{g.parts[1]}"])
        
        if gff.empty:
                gff = df
        else:
                gff = gff.merge(df, how = 'outer')
    gff = gff[gff['cond'].isin(['CDS', 'rRNA'])]
    gff = gff.T
    gff.columns = gff.iloc[0]
    gff = gff.iloc[1:]  

    assembly = assembly_stats.merge(gff, left_on = ['Name'], right_on= gff.index)
    cut = assembly['# Contigs'].mean() + (2* assembly['# Contigs'].std())
    assembly = assembly.rename(columns={'Name':'Isolate'})
    assembly['Quality'] = numpy.where(assembly['# Contigs'] <= cut, assembly['Quality'], 'WARNING - assembly outlier')
    assembly = assembly[['Isolate', 'bp','# Contigs','Ns','# Gaps','Min Contig size','Max Contig size','Avg Contig size','N50', 'CDS','rRNA', 'Quality' ]]

    assembly.to_csv(f"assembly.tab", sep = '\t', index = False)

    return assembly.to_dict(orient = 'records')

def open_toml(tml):

    data = toml.load(tml)

    return data

def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)
    
def main(assembly_stats, prokka):
    
    assemblies = assembly_combine(assemblies = assemblies)
    assembly_stat_dict = combine(prokka = prokka, assembly_stats = assemblies) 
    # set up data dict
    data = {}
    data['assembly_stat'] = assembly_stat_dict
    write_toml(data = data, output= f'assembly.toml')



if __name__ == '__main__':
    
    main(assembly_stats = f"{sys.argv[1]}", prokka = f"{sys.argv[2]}")
    

