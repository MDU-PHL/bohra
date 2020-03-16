import toml, pathlib, subprocess, sys, pandas, datetime
from snakemake import shell

def get_isolate_list(inputs):
    
    isolates = []
    for i in inputs:
        data = open_toml(i)
        isolate = pathlib.Path(i).parts[0]
        if data[isolate]['qc_snippy']['Quality'] == 'PASS':
            isolates.append(isolate)
    
    return ' '.join(isolates)

def core_stats(inputs):
    df = pandas.DataFrame()
    for i in inputs:
        tml = open_toml(i)
        isolate = list(tml.keys())[0]
        d = pandas.DataFrame(tml[isolate]['qc_snippy'], index = [0])
        d['Isolate'] = isolate
        d = d[['Isolate', 'Quality']]
        if df.empty:
            df = d
        else:
            df = df.append(d)

    # for core.txt
    tab = pandas.read_csv("core.txt", sep = '\t')
    tab['% USED'] = 100 * (tab['LENGTH'] - tab['UNALIGNED'])/ tab['LENGTH']
    tab['% USED'] = tab['% USED'].round(2)
    tab = tab.rename(columns={'ID':'Isolate'})
    tab = pandas.merge(left = df, right = tab, how = 'outer')
    tab = tab[['Isolate','LENGTH', 'ALIGNED','UNALIGNED','VARIANT','HET','MASKED','LOWCOV','% USED', 'Quality']]
    tab.to_csv(f"core_genome.tab", sep='\t', index = False)

    return df.to_dict(orient = 'records')


def generate_snippy_core_cmd(isolates, reference, mask):
    
    cmd = f"snippy-core --mask {mask} --ref {reference} {isolates}" if mask != 'nomask' else f"snippy-core --ref {reference} {isolates}"
    return cmd

def run_cmd(cmd):

    p = subprocess.run(cmd, shell = True, capture_output=True, encoding = 'utf-8')
    return p.returncode

def open_toml(tml):

    data = toml.load(tml)

    return data

def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)
    
def main(inputs, mask, reference):

    print('Getting a list of isolates to include in core analysis.')
    isolates = get_isolate_list(inputs)
    cmd = generate_snippy_core_cmd(isolates = isolates, reference = reference, mask = mask)
    print(f'Running core analysis with {cmd}')
    p = run_cmd(cmd)
    data = {}
    data['snippy-core'] = {}
    if p == 0:
        print('Core analysis has been performed and toml file will be updated.')
        date = datetime.datetime.today().strftime("%d_%m_%y")
        data['snippy-core'][date] = {}
        data['snippy-core'][date]['isolates'] = isolates.split()
        data['snippy-core'][date]['reference'] = reference
        data['snippy-core'][date]['mask'] = mask 
        data['snippy-core'][date]['vcf'] = 'core.vcf'
        data['snippy-core'][date]['aln'] = 'core.aln'
        data['snippy-core'][date]['txt'] = 'core.txt'
        data['snippy-core'][date]['tab'] = 'core.tab'
        data['snippy-core'][date]['fullaln'] = 'core.full.aln'
        data['snippy-core'][date]['stats'] = core_stats(inputs)

        write_toml(data = data, output = "snippy_core.toml")


#  {params.mask_string} {params.reference} {input}
mask_string = snakemake.params.mask_string
reference = snakemake.params.reference
inputs = snakemake.input

main(inputs = inputs, mask = mask_string, reference =reference)
    



# mash triangle -C *.msh

# mash sketch -m 5 -s 10000 -r -o 2019-12803-6/sketch -I 2019-12803-6 -C 2019-12803-6/R1.fq.gz 2019-12803-6/R1.fq.gz