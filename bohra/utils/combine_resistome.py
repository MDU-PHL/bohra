import toml, pathlib, subprocess, sys, pandas


def combine_resistome(inputs):

    tab = pandas.DataFrame()

    for i in inputs:
        tml = open_toml(i)
        isolate = tml.keys()[0]
        df = pandas.DataFrame(d = tml[isolate]['resistome']['data'], index = [0])
        df['Quality'] = 'PASS'
        if not tml[isolate]['resistome']['done']:
            df['Quality '] = 'NOT INCLUDED - FAIL QC'
        if tab.empty:
            tab = df
        else:
            tab = tab.append(df, sort = True)
    tab = tab.set_index('Isolate')
    tab.to_csv('resistome.tab', sep= '\t', index = True)
    return tab

def open_toml(tml):

    data = toml.load(tml)

    return data

def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)
    
def main(inputs):
    
    # set up data dict
    tab = combine_resistome(inputs)
    data = {}
    data['resistome'] = tab.to_dict(orient= 'records')
    write_toml(data = data, output= f'resistome.toml')



if __name__ == '__main__':
    
    main(inputs = f"{sys.argv[1]}")
    

