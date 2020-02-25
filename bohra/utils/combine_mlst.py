import toml, pathlib, subprocess, sys, pandas


def combine_mlst(inputs):

    tab = pandas.DataFrame()

    for i in inputs:
        tml = open_toml(i)
        isolate = tml.keys()[0]
        d = {
            "Isolate": isolate,
            "Scheme" : '',
            "ST": '',
            "Quality": 'PASS'
        }
        if tml[isolate]['mlst']['done'] and tml[isolate]['mlst']['data']['scheme'] != '-':
            for a in len(tml[isolate]['mlst']['data']['alleles']):
                d[f"Allele {a+1}"] = tml[isolate]['mlst']['data']['alleles'][a]
            d['Scheme'] = tml[isolate]['mlst']['data']['scheme']
            d['ST'] = tml[isolate]['mlst']['data']['ST']
            
        elif tml[isolate]['mlst']['done'] and tml[isolate]['mlst']['data']['scheme'] == '-':
            d["Scheme"] = "No scheme available"
        elif not tml[isolate]['mlst']['done']:
            df['Quality '] = 'NOT INCLUDED - FAIL QC'
        df = pandas.DataFrame(d, index = [0])
        if tab.empty:
            tab = df
        else:
            tab = tab.append(df, sort = True)
    tab = tab.set_index('Isolate')
    tab.to_csv('mlst.tab', sep= '\t', index = True)
    return tab

def open_toml(tml):

    data = toml.load(tml)

    return data

def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)
    
def main(inputs):
    
    # set up data dict
    tab = combine_mlst(inputs)
    data = {}
    data['mlst'] = tab.to_dict(orient= 'records')
    write_toml(data = data, output= f'mlst.toml')



if __name__ == '__main__':
    
    main(inputs = f"{sys.argv[1]}")
    

