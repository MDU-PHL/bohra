#!/usr/bin/env python3
import pathlib, subprocess, sys, datetime, pandas, re, numpy, jinja2, json


def write_tables(table, wd, job_id):
    '''
    Write a table, given a tab delimited file generate a html string
    '''
    # TODO add class isolate id to <tr>
    # TODO add class distances-isolate to tr if matrix and head-isolate to head td
    path =  f"{pathlib.Path(wd, job_id,table)}"
    
    data = open(path).readlines()
    # for header
    header = data[0].split('\t')
    
    if 'distances.tab' == table:
        tablehead = [f"<th class='{column}-head'>{column}</th>" for column in header]
    else:
        tablehead = [f"<th>{column}</th>" for column in header]
    # for body seqtablebody
    body = []
    for i in range(1,len(data)):
        raw = data[i].split('\t')
        raw = [r.strip() for r in raw]
        print(raw)
        if 'summary_table.tab' in table:
            row = [f"<tr class='{raw[0]} tiplab'>"]
        elif 'distances.tab' in table:
            row = [f"<tr class='distances-{raw[0]}'>"]
        elif 'assembly' in table:
            row = [f"<tr class='{raw[0]}-assembly'>"]
        elif 'species' in table:
            row = [f"<tr class='{raw[0]}-species-identification'>"]
        elif 'core_genome' in table:
            row = [f"<tr class='{raw[0]}-core-genome'>"]
        elif 'mlst' in table:
            row = [f"<tr class='{raw[0]}-mlst'>"]
        elif 'resistome' in table:
            row = [f"<tr class='{raw[0]}-resistome'>"]
        elif 'virulome' in table:
            row = [f"<tr class='{raw[0]}-virulome'>"]
        elif 'seqdata' in table:
            row = [f"<tr class='{raw[0]}-sequence-data'>"]
        else:
            row = [f"<tr>"] # TODO add class isolate id to <tr>
        if 'distances.tab' == table:
                for d in range(len(raw)):
                    row.append(f"<td align=\"center\" class = \"{raw[0]}_{header[d].strip()}\">{raw[d]}</td>")    
        elif 'resistome' in table:
            for d in raw:
                # 
                dr = d.split(';')
                drs = []
                for i in dr:
                    if '^' in i:
                        drs.append(f"<span style=\"color:#3483eb\">{i.strip('^')}</span>")
                    else:
                        drs.append(i)
                x = ';'.join(drs)
                row.append(f"<td align=\"center\">{x}</td>")
        else:
            for d in raw:
                row.append(f"<td align=\"center\">{d}</td>")
        row.append(f"</tr>")
        body = body + row
        
    return('\n'.join(tablehead),'\n'.join(body))

def get_isolates_preview(table):
    
    data = open(table).readlines()
    isolates = []
    for d in range(1,len(data)):
        l = data[d].split()
        isolates.append(l[0])
    return len(isolates)

def preview_distances_tab(table):
    
    data = open(table).readlines()
    length = get_isolates_preview(table)
    header = ['mash-dist']
    
    body = []
    for d in range(1,len(data)):
        row = [f"<tr>"]
        l = data[d].split()
        header.append(l[0])
        for i in range(len(l)):
        # print(d)  
            row.append(f"<td align=\"center\" class = \"{l[0]}_{header[i]}\">{l[i]}</td>")
        row.append(f"</tr>")
        body = body + row
    # print(body)
    tablehead = [f"<th class='{column}-head'>{column}</th>" for column in header]
    return('\n'.join(tablehead),'\n'.join(body))

def get_tree_string(pipeline, wd, job_id):
    '''
    Generate a tree image from a newick
    input:
        :reportdir: the directory where report files are stored
    output:
        string reporesentation of the path to the tree image
    '''
    tree_file = f"{pathlib.Path(wd, job_id,'preview.newick')}" if pipeline == 'preview' else f"{pathlib.Path(wd, job_id,'core.treefile')}"
    with open(f"{tree_file}", 'r') as t:
        tree = t.read().strip()

    return tree
def adjust_offset(row, d):
    if row['CHR'] in d:
        return(int(d[row['CHR']]) + int(row['POS']))
    else:
        return(int(row['POS']))

def generate_dict(idx_file):

    d = {}
    pos = 0
    offset = 0
    with open(f"{idx_file}") as f:
        for line in f:
            l = line.split()
            d[l[0]] = offset
            offset = offset + int(l[1])
    return(d)

def plot_snpdensity():

    '''
    generate a snp-density accross the genome plot - using the core.tab file
    input:
        :reportdir: the directory where files are kept
    out:
        :snpdensityscript: the javascript string for insert into html
        :snpdenstiydiv: the html div for inster in html doc
    '''
    # helper functions for getting the data into the right format.
    
    
    # open fai file and generate the dictionary
    idx = pathlib.Path('ref.fa.fai')
    d = generate_dict(idx)
    # get the core file
    core = 'core.tab'
    df = pandas.read_csv(core, sep = '\t')
    # get a list of isolate names
    names = list(df.columns[3:len(df.columns)])
    # if the there is no snp in the isolate (ie same as ref then mak na - then easy to drop)
    for i in names:
        df[i]=numpy.where(df['REF'] == df[i], numpy.nan, df[i])
    # generate the offset value
    df['POS_OFFSET'] = df[['CHR', 'POS']].apply(lambda x:adjust_offset(x,d), axis = 1)
    # collect positions to get allow for histogram and dropna (no snp)
    melted_df = pandas.melt(df, id_vars=['POS_OFFSET'], value_vars=names)
    melted_df = melted_df.dropna()
    melted_df = melted_df.sort_values(by= ['POS_OFFSET'])
    # generate histogram
    # snpdensityscript, snpdensitydiv = self.plot_histogram(series=melted_df['POS']/1000,xlabel="Genome Position (MB)", ylabel="SNPS",bins=10000)
    contig_breaks = [d[value] for value in d]
    
    # return dictionary
    return(list(melted_df['POS_OFFSET']/1000))


def plot_distances():

    '''
    generate a snp-density plot - using the distacnes.tab file
    input:
        :reportdir: the directory where files are kept
    out:
        :distancescript: the javascript string for insert into html
        :distancesdiv: the html div for inster in html doc
    '''
    distance = 'distances.tab'

    df = pandas.read_csv(distance, sep = '\t')
    # get a list of isolate names
    names = list(df.columns[1:len(df.columns)])
    col1 = df.columns[0]
    
    # if the there is no snp in the isolate (ie same as ref then mak na - then easy to drop)
    
    # collect positions to get allow for histogram and dropna (no snp)
    melted_df = pandas.melt(df, id_vars=[col1], value_vars=names)
    melted_df = melted_df[melted_df[col1]!= melted_df['variable']]
    # generate histogram
    # distancescript, distancediv = self.plot_histogram(series=melted_df['value'], xlabel="SNP distances", ylabel="Frequency", bins=100)
    # td['pairwisedistance'] = {'script':distancecript, 'div':distancediv}
    # return dictionary
    return(list(melted_df['value']))

def get_preview_dict(run_kraken):
    print(run_kraken)
    td = [{'title': 'Preview', 'link':'preview-tree', 'file': 'preview.newick', 'type': 'tree'},{'file':'seqdata.txt', 'title':'Sequence Data', 'link': 'sequence-data', 'type' : 'table'}]
    if run_kraken == 'speciation':
        td.append({'file':'species.txt', 'title': 'Species Identification', 'type': 'table', 'link':'species-identification'})
    return td

def snps_dict(td):
    
    core_genome_td = {'file': 'core_genome.tab', 'title': 'Core Genome', 'link':'core-genome', 'type':'table'}
    snp_density_td = {'title': 'SNP density', 'link':'snp-density', 'type':'graph'}
    core_phylogeny_td = {'title': 'Phylogeny', 'link':'phylogeny', 'file': 'core.treefile', 'type': 'tree'}
    snp_distance_td = {'file': 'distances.tab', 'title':'SNP distances', 'type':'matrix', 'link':'snp-distances'}
    snps = [core_genome_td, snp_density_td, core_phylogeny_td,snp_distance_td]
    for s in snps:
        td.append(s)
    
    return td

def assembly_dict(td):
    for_td = []
    mlst_td = {'file':'mlst.txt', 'title':'MLST', 'type':'table', 'link':'mlst'}
    resistome_td = {'file':'resistome.txt', 'title':'Resistome', 'type':'table', 'link':'resistome'}
    # list of assembly tasks
    assembly_stat_td = {'file': 'assembly.txt', 'title':'Assembly', 'type':'table', 'link':'assembly'}
    a_td = [assembly_stat_td, mlst_td, resistome_td]
    for a in a_td:
        if pathlib.Path(a['file']).exists():
            for_td.append(a)
    for f in for_td:
        td.append(f)

    return td

def sa_dict(td):

    td = snps_dict(td)
    td = assembly_dict(td)

    return td

def roary_dict(td):

    roary_td = {'file':'roary/summary_statistics.txt', 'title':'Pan Genome', 'type': 'pan', 'image': f"{pathlib.Path('pan_genome.svg').open().read()}", 'link':'pan-genome'}

    td.append(roary_td)

    return td

def all_dict(td):

    td = sa_dict(td)
    td = roary_dict(td)

    return td

def get_dict(pipeline, run_kraken):

    td = [{'file':'seqdata.txt', 'title':'Sequence Data', 'link': 'sequence-data', 'type' : 'table'}, {'file':'summary_table.txt','title':'Summary', 'link':'summary', 'type':'summary'},{'file': 'software_versions.txt', 'title': 'Tools', 'type': 'versions', 'link':'versions'}]
    
    if run_kraken == 'speciation':
        td.append({'file':'species.txt', 'title': 'Species Identification', 'type': 'table', 'link':'species-identification'})
    if pipeline == 's':
        td = snps_dict(td)
    elif pipeline == 'a':
        td = assembly_dict(td)
    elif pipeline == 'sa':
        td = sa_dict(td)
    elif pipeline == 'all':
        td = all_dict(td)
    # print(td)
    return td

def fill_vals(td, pipeline, wd, job_id):

    for t in range(len(td)):
    
        # TODO if table add a modal modal + link and link will be title lowercase with hyphen
        if td[t]['type'] == 'table':
            td[t]['head'], td[t]['body'] = write_tables(table=td[t]['file'], wd = wd, job_id=job_id)
        if td[t]['type'] == 'pan':
            td[t]['head'], td[t]['body'] = write_tables(table=td[t]['file'], wd = wd, job_id=job_id)
        if td[t]['type'] == 'matrix':
            td[t]['head'], td[t]['body'] = write_tables(table=td[t]['file'], wd = wd, job_id=job_id)
            # snpdistances = plot_distances()
        # if td[t]['link'] == 'snp-density':
            # snpdensity= plot_snpdensity()
        if td[t]['type'] == 'versions':
            td[t]['head'], td[t]['body'] = write_tables(table=td[t]['file'], wd = wd, job_id=job_id)
        if td[t]['type'] == 'summary':
            # generate_summary()
            td[t]['head'], td[t]['body'] = write_tables(table = td[t]['file'], wd = wd, job_id=job_id)
    return td

def get_software_versions(software):

    '''
    Given the name of the software, find the version
    input:
        :software: the name of the software
    output:
        a string in the form of 'Name_of_Sofware v.X.Y.Z'
    '''
    
    version_pat = re.compile(r'\bv?(?P<major>[0-9]+)\.(?P<minor>[0-9]+)\.(?P<release>[0-9]+)(?:\.(?P<build>[0-9]+))?\b')

    if software == 'snp-dists':
        vs = '-v'
    else:
        vs = '--version'
    cmd = f"{software} {vs} 2>&1"
    print(f"Checking software with {cmd}")
    p = subprocess.run(cmd, shell = True, capture_output=True, encoding = "utf-8")
    sft = p.stdout
    
    v = version_pat.search(sft)
    v = v.group()
    sft_version = f"{software} v.{v}"
    return(sft_version)

def make_dict_versions(tools):
    '''
    Called by get_software_file to make a dictionary of tools used
    input:
        :tools: a list of tools
    output:
        a dictionary with tools and versions.
    '''
    tool_dict = {}
    for t in tools:
            v = get_software_versions(t)
            tool_dict[t] = v
    return(tool_dict)


def get_software_file(pipeline, assembler = ''):
    '''
    get the versions of software on the system at completion of the piepline
    input:
        :reportdir: the directory where report files are stored
        :pipeline: the type of pipeline
        :assembler: the assembler used in the pipeline
    '''
    
    
    p = pathlib.Path('software_versions.tab')

    p.write_text('\n'.join(versions))

def merge_dfs(start, added):
    if start.empty:
        start = added
    else:
        start = start.merge(added, how = 'outer')
    return(start)

def generate_summary(wd, job_id):
    '''
    function to generate a summary table
    '''
    p = pathlib.Path(wd, job_id, 'report')
    tabs = [t for t in p.iterdir() if f"{t.suffix}" == '.txt']
    print(tabs)
    summary_df = pandas.DataFrame()
    df_list = []
    
    # print(tabs)
    for tab in tabs:
        # print(df)
        print(tab)
        if 'species' in f"{tab.name}":
            species = pandas.read_csv(tab, sep = '\t')
            print(species)
            species = species[['Isolate', 'Match 1']]
            summary_df = merge_dfs(summary_df, species)
            summary_df = summary_df.rename(columns={'Match 1': 'Species'})
        elif 'seqdata' in f"{tab.name}":
            seq = pandas.read_csv(tab, sep = '\t')
            print(seq)
            seq = seq[['Isolate', 'Estimated average depth']]
            
            summary_df = merge_dfs(summary_df, seq)
        elif 'assembly' in f"{tab.name}":
            assembly = pandas.read_csv(tab, sep = '\t')
            assembly = assembly[['Isolate', '# Contigs']]
            print(assembly)
            summary_df = merge_dfs(summary_df, assembly)
        elif 'mlst' in f"{tab.name}":
            mlst = pandas.read_csv(tab, sep = '\t', header=None)
            print(mlst)
            mlst = mlst.rename(columns = {0:'Isolate', 2:'ST'})
            mlst = mlst[['Isolate', 'ST']]
            print(mlst)
            summary_df = merge_dfs(summary_df, mlst)
        elif 'core_genome' in f"{tab.name}":
            core = pandas.read_csv(tab, sep = '\t')
            print(core)
            core = core[['Isolate', '% USED']]
            summary_df = merge_dfs(summary_df, core)
    print(summary_df)
    isolates = len(summary_df['Isolate'])
    # print(mlst)
    # print(summary_df)
    summary_df = summary_df.fillna('NA')
    summary_file = 'summary_table.txt'
    summary_df.to_csv(summary_file, sep = '\t', index = False)

    return isolates

def return_tables(pipeline):

    if pipeline == 'preview':
        
        tables =['sequence-data','species-identification']
        modaltables =['sequence-data','species-identification']
        display = f""
    
    elif pipeline == 'default':
        
        tables =['core-genome', 'snp-distances', 'mlst', 'assembly', 'resistome', 'sequence-data','species-identification']
        modaltables = ['core-genome',  'mlst', 'assembly', 'resistome', 'sequence-data','species-identification']
        display = f""
        # td.extend(s_td)
    elif pipeline == 'all':
        # a_ll = td.extend()
        tables =['core-genome', 'species-identification','snp-distances', 'mlst', 'assembly', 'resistome', 'sequence-data', 'pan-genome']
        modaltables = ['core-genome',  'mlst', 'assembly', 'resistome', 'sequence-data']
        display = f""
    
    tables.append('versions')

    return tables, modaltables, display

    
def main(wd, pipeline,job_id, resources, run_kraken, day, user,assembler = ''):

    # get analysis dict
    _dict = json.load(open(f"{pathlib.Path(resources,'report_analysis.json')}"), 'r')
    
    p = pathlib.Path('.')
    print(p)
    reporthtml = pathlib.Path('report.html')
    print(reporthtml)
    # path to html template
    indexhtml = pathlib.Path(resources,'index.html') # replace with template
    print(indexhtml)
    # initialise dictionary
    print(pipeline)
    # print(inputs)
    data = {
        'newick' :'',
        'snpdensity':'',
        'snpdistances':'',
        'display':'',
        'job_id':job_id,
        'pipeline':pipeline,
        'date':day,
        'user':user,
        'tree_heigth':0,
        'modaltables':'',
        'tables':''
        }
    # print(data)
    if pipeline == 'preview':
        td = get_preview_dict(run_kraken = run_kraken)
        data['tree_height'] = get_isolates_preview(f"{pathlib.Path(wd,job_id,'preview_distances.tab')}") * 25
        tables, modaltables, display = return_tables(pipeline = pipeline)
        data['tables'] = tables
        data['modaltables'] = modaltables
        data['display'] = display
        # print(td)
    else:
        isos = generate_summary(wd = wd, job_id = job_id)
        # print(isos)
        # get_software_file(pipeline = pipeline, assembler = assembler)  
        td = get_dict(pipeline = pipeline, run_kraken = run_kraken)
        tables, modaltables, display = return_tables(pipeline = pipeline)
        data['tables'] = tables
        data['modaltables'] = modaltables
        data['display'] = display
        data['tree_heigth'] = isos * 25

    if pipeline not in ['preview']:
        data['snpdensity']= plot_snpdensity()
        data['snpdistances']= plot_distances()
        data['newick'] = get_tree_string(pipeline = pipeline, wd = wd, job_id = job_id)
    elif pipeline in ['default', 'preview', 'all']:
        data['newick'] = get_tree_string(pipeline = pipeline, wd = wd, job_id= job_id)

    

# newick = newick, display = display,tables = tables,td = td, job_id = job_id, pipeline = pipeline, snpdistances=snpdistances, snpdensity = snpdensity, modaltables = modaltables, date = date
    td = fill_vals(td=td, pipeline = pipeline, wd = wd, job_id=job_id)
    data['td'] = td
    print("rendering html")
    report_template = jinja2.Template(pathlib.Path(indexhtml).read_text())
    reporthtml.write_text(report_template.render(data))
    

pipeline = sys.argv[1]
job_id = sys.argv[2]
resources = sys.argv[3]
wd = sys.argv[4]
day = sys.argv[6]
user = sys.argv[7]
run_kraken = 'speciation' if sys.argv[5]=='true' else ''

if pipeline != 'preview':
    assembler = sys.argv[5]
else:
    assembler = ''
main(pipeline = pipeline,wd = wd, job_id= job_id, 
    assembler= assembler, day = day, user=user, 
    resources = resources, run_kraken = run_kraken)
