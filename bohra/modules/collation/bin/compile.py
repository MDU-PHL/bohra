#!/usr/bin/env python3
import pathlib, subprocess, sys, datetime, pandas, re, numpy, jinja2, json, csv
from Bio import SeqIO

def _write_tables(table, wd, job_id, link):
    '''
    Write a table, given a tab delimited file generate a html string
    '''
    # TODO add class isolate id to <tr>
    # TODO add class distances-isolate to tr if matrix and head-isolate to head td
    path =  f"{pathlib.Path(wd, job_id,'report', table)}"
    
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
        if 'distances' in table:
            row = [f"<tr class='distances-{raw[0]}'>"]
        else:
            row = [f"<tr class='{raw[0]}-{link}'>"]
        if 'distances.tab' == table:
            # to allow for display on tree
                for d in range(len(raw)):
                    row.append(f"<td align=\"center\" class = \"{raw[0]}_{header[d].strip()}\">{raw[d]}</td>")    
        elif 'resistome' in table:
            # to allow for coloring of partials
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
        else: #default
            for d in raw:
                row.append(f"<td align=\"center\">{d}</td>")
        row.append(f"</tr>")
        body = body + row
        
    return('\n'.join(tablehead),'\n'.join(body))

def _get_tree_string(pipeline, wd, job_id):
    '''
    Generate a tree image from a newick
    input:
        :reportdir: the directory where report files are stored
    output:
        string reporesentation of the path to the tree image
    '''
    tree_file = pathlib.Path(wd, job_id,'preview.newick') if pipeline == 'preview' else pathlib.Path(wd, job_id,'core.treefile')
    if tree_file.exists():
        with open(f"{tree_file}", 'r') as t:
            tree = t.read().strip()
    else:
        tree = 'No tree available'
    return tree


def _get_offset(reference):

    d = {}
    offset = 0
    records = list(SeqIO.parse(reference, "genbank"))
    if records == []:
        records = list(SeqIO.parse(reference, "fasta"))
    if records != []:
        for record in records:
            d[record.id] = {'offset' : offset, 'length': len(record.seq)}
            offset += len(record.seq)
    
    return d, offset

def _plot_snpdensity(reference,wd, job_id, isos):

    '''
    generate a snp-density accross the genome plot - using the core.tab file
    input:
        :reportdir: the directory where files are kept
    out:
        :snpdensityscript: the javascript string for insert into html
        :snpdenstiydiv: the html div for inster in html doc
    '''
    # helper functions for getting the data into the right format.
    
    _dict,offset = _get_offset(reference = reference)
    _all_pos = list(range(1,offset+1))
    _snp_dict = {}
    # collate all snps in snps.tab
    for i in isos:
        # open snps.tab
        snps = pathlib.Path(wd, job_id, 'report', i, 'snps.tab')
        if snps.exists():
            with open(snps, 'r') as s:
                reader = csv.DictReader(s, delimiter = '\t')
                for row in reader:
                    chrom = row['CHROM']
                    ofs = _dict[chrom]['offset'] # offset value of this chromosome
                    pos = int(row['POS']) + ofs #get the position in the genome (with offset)
                    if pos in _snp_dict: 
                        # if the pos is in the dict already, it has been found in another sample so increment
                        _snp_dict[pos] = _snp_dict[pos] + 1
                    else:
                        _snp_dict[pos] = 1
    # now generate list for x value in graph
    _density = []
    for a in _all_pos:
        if a in _snp_dict:
            _density[a].append(_snp_dict[a])
        else:
            _density[a].append(0)
    # open fai file and generate the dictionary
    
    
    # return dictionary
    return(_density)


def _plot_distances():

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
    
    return(list(melted_df['value']))

def _get_pan_genome(image, wd, job_id):
    path = pathlib.Path(wd, job_id, 'report', image)
    if path.exists():
        with open(f"{path}", 'r') as f:
            return f.read().strip()
    else:
        return ''

def _fill_vals(td, pipeline, wd, job_id):

    for t in range(len(td)):
    
        # TODO if table add a modal modal + link and link will be title lowercase with hyphen
        if td[t]['type'] == 'table':
            td[t]['head'], td[t]['body'] = _write_tables(table=td[t]['file'], wd = wd, job_id=job_id, link = td[t]['link'])
        if td[t]['type'] == 'pan':
            td[t]['head'], td[t]['body'] = _write_tables(table=td[t]['file'], wd = wd, job_id=job_id, link = td[t]['link'])
            td[t]['image'] = _get_pan_genome(image = td[t]['file'], wd = wd, job_id=job_id)
        if td[t]['type'] == 'matrix':
            td[t]['head'], td[t]['body'] = _write_tables(table=td[t]['file'], wd = wd, job_id=job_id, link = td[t]['link'])
            # snpdistances = plot_distances()
    return td



def _return_tables(pipeline):

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

def _get_isos(isolates):

    with open(isolates, 'r') as f:
        isos = f.read().strip().split('\n')
    
    return isos

def main(args):

    # get analysis dict
    _dict = json.load(open(f"{pathlib.Path(args.template_dir, 'report_analysis.json')}", 'r'))
    isos = _get_isos(args.isolates)
    reporthtml = pathlib.Path('report.html')
    # path to html template
    indexhtml = pathlib.Path(args.template_dir,'index.html') 
    # initialise dictionary
    # data is the dictionary passed to jinja2 to fill html
    data = {
        'newick' :'',
        'snpdensity':'',
        'snpdistances':'',
        'display':'',
        'job_id':args.job_id,
        'pipeline':args.pipeline,
        'date':args.day,
        'user':args.user,
        'tree_heigth':0,
        'modaltables':'',
        'tables':''
        }
    # print(data)
    td = _dict[pipeline]
    tables, modaltables, display = _return_tables(pipeline = args.pipeline)
    data['tree_height'] = len(isos) * 25
    data['tables'] = tables
    data['modaltables'] = modaltables
    data['display'] = display
    data['newick'] = _get_tree_string(pipeline = args.pipeline, wd = args.launchdir, job_id = args.job_id)
        # print(td)
    # generate_summary(wd = wd, job_id = job_id)
        # print(isos)
        # get_software_file(pipeline = pipeline, assembler = assembler)  
    if pipeline not in ['preview']:
        data['snpdensity']= _plot_snpdensity(reference = args.reference,wd = args.launchdir, job_id = args.job_id, isos = isos)
        
    if isos > 1:
        data['snpdistances']= _plot_distances()
    

# newick = newick, display = display,tables = tables,td = td, job_id = job_id, pipeline = pipeline, snpdistances=snpdistances, snpdensity = snpdensity, modaltables = modaltables, date = date
    td = _fill_vals(td=td, pipeline = pipeline, wd = wd, job_id=job_id)
    data['td'] = td
    print("rendering html")
    report_template = jinja2.Template(pathlib.Path(indexhtml).read_text())
    reporthtml.write_text(report_template.render(data))
    

def set_parsers():
    # setup the parser
    parser = argparse.ArgumentParser(description='Collate and write bohra report',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--pipeline',
        help='',
        default = 'default')
    parser.add_argument('--launchdir',
        help='',
        default = '')
    parser.add_argument('--template_dir',
        help='',
        default = '')
    parser.add_argument('--day',
        help=f'',
        default = '')
    parser.add_argument('--user', 
        help='',
        default = '')
    parser.add_argument('--job_id',
        help='',
        default = '')
    parser.add_argument('--isolates',
        help = '',
        default = '')
    parser.add_argument('--reference',
        help = '',
        default = '')
    
    
    
    parser.set_defaults(func=main)
    args = parser.parse_args()
    
    if vars(args) == {}:
        parser.print_help(sys.stderr)
    else:
        # print(args)
        args.func(args)
	

def main():
    """
    run pipeline
    """

    args = set_parsers()
    

if __name__ == "__main__":
    main()

