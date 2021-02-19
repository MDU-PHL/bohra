import toml, pathlib, subprocess, sys, pandas, os

from snakemake import shell

def get_top_3(isolate):
    
    report = pathlib.Path(isolate, 'kraken2.tab')
    df= pandas.read_csv(report, sep = "\t", header =None, names = ['percentage', 'frag1', 'frag2','code','taxon','name'])
    df['percentage'] = df['percentage'].apply(lambda x:float(x.strip('%')) if isinstance(x, str) == True else float(x)) #remove % from columns
    df = df.sort_values(by = ['percentage'], ascending = False)
    df = df[df['code'].isin(['U','S'])]     
    df = df.reset_index(drop = True) 
    top_three = [
        (f"{df.loc[0,'name'].strip()}",round(df.loc[0,'percentage'],2)),
        (f"{df.loc[1,'name'].strip()}",round(df.loc[1,'percentage'],2)),
        (f"{df.loc[2,'name'].strip()}",round(df.loc[2,'percentage'],2))
    ]

    return top_three

def generate_cmd(prefill, r1, r2, isolate, db, threads, data):
    
    prfl = pathlib.Path(prefill, isolate, 'kraken2.tab')
    if prfl.exists() and os.access(f"{prfl}", os.R_OK):
        cmd = f"cp {prfl} {isolate}/kraken2.tab"
        data[isolate]['kraken']['kraken_db'] = f"Prefilled from: {prfl}"
        print(f"Species id will be determined from results retrieved from {prfl}")
    else:
        cmd = f"kraken2 --paired {r1} {r2} --threads {threads} --minimum-base-quality 13 --report {isolate}/kraken2.tab --memory-mapping --db {db}"
        data[isolate]['kraken']['kraken_db'] = f"{db}"
        print(f"Species id will be determined from {db}")
    print(cmd)
    return cmd, data

def run_cmd(cmd):
    
    p = subprocess.run(cmd, shell = True, capture_output=True, encoding = 'utf-8')
    return p.returncode
    

def open_toml(tml):

    data = toml.load(tml)

    return data

def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)
    
def main(r1, r2, isolate, kraken_db,prefill, kraken_threads):
    
    # set up data dict
    tml = f"{pathlib.Path(isolate, 'kraken.toml')}"
    data = {}
    data[isolate] = {}
    data[isolate]['kraken'] = {}
    # run kraken
    print(f"Identifying species of isolate {isolate}")
    cmd, data = generate_cmd(prefill = prefill, r1 = r1, r2 = r2, isolate = isolate, db = kraken_db, threads = kraken_threads, data = data)
    kraken_returncode = run_cmd(cmd)
    if kraken_returncode == 0:
        # add to data dict
        print(f"Kraken successfully ran - you will be able to see the top 3 species id in the toml file.")
        data[isolate]['kraken']['done'] = True
        data[isolate]['kraken']['file'] = f"{isolate}/kraken2.tab"
        data[isolate]['kraken']['top_3'] = get_top_3(isolate)
        
        write_toml(data = data, output= f"{isolate}/kraken.toml")
    else:
        print(f"Somthing went wrong with species identification, please check your inputs.")


print("Running kraken")

r1 = snakemake.input.r1
# r1 = 'test2/2012-09771-B/R1.fq.gz'
r2 = snakemake.input.r2
# r2 = 'test2/2012-09771-B/R2.fq.gz'
isolate = snakemake.wildcards.sample
# isolate = '2012-09771-B'
kraken_db = snakemake.params.kraken_db
# kraken_db = '/home/linuxbrew/db/kraken2/microbe'
kraken_threads = snakemake.threads
# kraken_threads = 16
prefill = snakemake.params.prefill_path
# prefill = ''
# print(kraken_threads)

main(r1 = r1, r2 = r2, isolate = isolate, kraken_db = kraken_db, kraken_threads = kraken_threads, prefill = prefill)


# if __name__ == '__main__':
    
#     main(r1 = f"{sys.argv[1]}", r2 = f"{sys.argv[2]}", isolate = f"{sys.argv[3]}", kraken_db = f"{sys.argv[4]}", prefill = f"{sys.argv[5]}")
    
# {input.r1} {input.r2} {wildcards.sample} {params.kraken_db} {params.prefill_path}
