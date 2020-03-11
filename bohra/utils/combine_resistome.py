import toml, pathlib, subprocess, sys, pandas
from snakemake import shell

def combine_resistome(inputs):

    tab = pandas.DataFrame()

    for i in inputs:
        tml = open_toml(i)
        isolate = list(tml.keys())[0]

        if tml[isolate]['resistome']['done']:
            df = pandas.DataFrame(tml[isolate]['resistome']['data'], index = [0])
            df['Quality'] = 'PASS'
        else:
            df = pandas.DataFrame({'Isolate':isolate, 'Quality': 'NOT INCLUDED - FAIL QC'}, index = [0])
            # df['Quality '] = 'NOT INCLUDED - FAIL QC'
        if tab.empty:
            tab = df
        else:
            tab = tab.append(df, sort = True)
    tab = tab.set_index('Isolate')
    test = tab[tab['Quality'] == 'PASS']
    if not test.empty:
        tab.to_csv('resistome.tab', sep= '\t', index = True)
        return tab
    else:
        return tab

def open_toml(tml):

    data = toml.load(tml)

    return data

def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)
    
def main(inputs):
    
    # set up data dict
    
    data = {}
    tab = combine_resistome(inputs)
    if tab.empty:
        data['resistome'] = 'resistome not performed'
    else:
        data['resistome'] = tab.to_dict(orient= 'records')
    write_toml(data = data, output= f'resistome.toml')

inputs = snakemake.input

main(inputs = inputs)
    

