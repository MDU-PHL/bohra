#!/usr/bin/env python3
import pathlib, subprocess, argparse, datetime, pandas, re, numpy, jinja2, json, csv
import altair as alt
from Bio import SeqIO

alt.data_transformers.disable_max_rows()

def _get_tree_string(pipeline, wd,phylo):
    '''
    Generate a tree image from a newick
    input:
        :reportdir: the directory where report files are stored
    output:
        string reporesentation of the path to the tree image
    '''
    tree_file = pathlib.Path(wd, 'report','preview.newick') if pipeline == 'preview' else pathlib.Path(wd, 'report','core.newick')
    if tree_file.exists() and phylo == 'true':
        with open(f"{tree_file}", 'r') as t:
            tree = t.read().strip()
    else:
        tree = 'No tree available'
    print(tree)
    return tree


def _get_offset(reference):
    # print(reference)
    d = {}
    offset = 0
    records = list(SeqIO.parse(reference, "genbank"))
    if records == []:
        records = list(SeqIO.parse(reference, "fasta"))
    if records != []:
        for record in records:
            d[record.id.split('.')[0]] = {'offset' : offset, 'length': len(record.seq)}
            offset += len(record.seq)
    # print(d)
    return d, offset

def _plot_snpdensity(reference,wd, isos):

    '''
    generate a snp-density accross the genome plot - using the core.tab file
    input:
        :reportdir: the directory where files are kept
    out:
        :snpdensityscript: the javascript string for insert into html
        :snpdenstiydiv: the html div for inster in html doc
    '''
    # helper functions for getting the data into the right format.
    
    _dict,offset = _get_offset(reference = f"{pathlib.Path(wd,reference)}")
    _all_pos = list(range(1,offset+1))
    # print(_dict)
    _snp_dict = {}
    # collate all snps in snps.tab
    for i in isos:
        # open snps.tab
        snps = pathlib.Path(wd, i, 'snps.tab')
        if snps.exists():
            with open(snps, 'r') as s:
                reader = csv.DictReader(s, delimiter = '\t')
                for row in reader:
                    # print(row)
                    chrom = row['CHROM'].split('.')[0]
                    ofs = _dict[chrom]['offset'] # offset value of this chromosome
                    pos = int(row['POS']) + ofs #get the position in the genome (with offset)
                    if pos in _snp_dict: 
                        # print(i)
                        # if the pos is in the dict already, it has been found in another sample so increment
                        _snp_dict[pos] = _snp_dict[pos] + 1
                    else:
                        # print(i)
                        _snp_dict[pos] = 0
    # now generate list for x value in graph
    _density = []
    
    for a in _all_pos:
        
        if a in _snp_dict:
            _density.append(_snp_dict[a])
        else:
            _density.append(0)
    # print(max(_snp_dict.values()))
    # print(_snp_dict)
    # return dictionary
    # _snp_dict = {1:5,3:3,10:2}
    _df = pandas.DataFrame.from_dict(_snp_dict, orient='index').reset_index()
    _df = _df.rename(columns={0:'snps', 'index':'Genome_position'})
    # _df = _df[_df['snps'] != 0]
    # print(_df)
    bins = list(range(1,max(_all_pos),1000))
    # print(list(b))
    s = _df.groupby(pandas.cut(_df['Genome_position'], bins=bins)).size()
    df  = s.to_frame().reset_index().reset_index()
    df['index'] = df['index'].apply(lambda x: (x + 50)*1000)
    df = df.rename(columns = {0:'snps'})
    df = df[['index','snps']]
    
    chart = alt.Chart(df).mark_bar().encode(
                            x=alt.X('index', axis=alt.Axis(title='Genome position'), scale = alt.Scale(domain=(1, max(bins)))),
                            y=alt.Y('snps', axis = alt.Axis(title = "SNPs (per 1000 bases)"))
                        ).properties(
                            width = 1200,
                            height = 200
                        )
    chart = chart.to_json()
    return chart

def _plot_distances(wd):

    '''
    generate a snp-density plot - using the distacnes.tab file
    input:
        :reportdir: the directory where files are kept
    out:
        :distancescript: the javascript string for insert into html
        :distancesdiv: the html div for inster in html doc
    '''
    distance = f"{pathlib.Path(wd, 'report','distances.tab')}"
    try:
        df = pandas.read_csv(distance, sep = '\t')
        # get a list of isolate names
        names = list(df.columns[1:len(df.columns)])
        col1 = df.columns[0]
        
        # if the there is no snp in the isolate (ie same as ref then mak na - then easy to drop)
        
        # collect positions to get allow for histogram and dropna (no snp)
        melted_df = pandas.melt(df, id_vars=[col1], value_vars=names)
        melted_df = melted_df[melted_df[col1]!= melted_df['variable']]
        # brush = alt.selection_interval(encodings=['x'])
        # base = alt.Chart(melted_df).mark_bar().encode(
        #                         y='count():Q'
        #                     ).properties(
        #                         width=1200,
        #                         height=100
        #                     )
        # chart = alt.vconcat(
        #             base.encode(
        #                 alt.X('value:Q',
        #                 bin=alt.Bin(maxbins=30, extent=brush),
        #                 scale=alt.Scale(domain=brush)
        #                 )
        #             ),
        #             base.encode(
        #                 alt.X('value:Q', bin=alt.Bin(maxbins=30)),
        #             ).add_selection(brush)
        #             )
        
        # chart = chart.to_json()
        chart = alt.Chart(melted_df).mark_bar().encode(
                            alt.X('value', axis = alt.Axis(title = 'Pairwise SNP distance')),
                            y=alt.Y('count()', axis= alt.Axis(title = "Frequency"))
                        ).properties(
                                width=1200,
                                height=200
                            )
        chart = chart.to_json()
        return chart
    except:
        return {}

def _get_pan_genome(d, wd):

    image = d[1]['image']
    path = pathlib.Path(wd, 'report', image)
    # print(path)
    if path.exists():
        with open(f"{path}", 'r') as f:
            return f.read().strip()
    else:
        return ''

def _get_isos(wd, iso_list):
    p = pathlib.Path(wd, iso_list)
    with open(p, 'r') as f:
        isos = f.read().strip().split('\n')
    
    return isos

def _get_versions(wd):

    p = pathlib.Path(wd, 'report', 'software_versions.txt')
    if p.exists():
        with open(f"{p}", 'r') as f:
            data = f.read().strip().split('\n')
            _head = f"<th class='version-head'>{data[0]}</th>"
            body = []
            for d in data[1:]:
                body.append(f"<tr><td>{d}</td></tr>")
        return _head,'\n'.join(body)
    else:
        return f"<th class='version-head'>Nothing to display</th>",""

def _generate_table(d, columns,comment, tables, wd, iso_dict, id_col):

    if id_col == '':
        return tables,columns,comment
    cols = []
    with open(f"{pathlib.Path(wd, 'report', d['file'])}", 'r') as f:
        reader = csv.DictReader(f, delimiter = '\t')
        comment[d['link']] = d['comment']
        tables[d['link']] = {'table':[], 'name': d['title'], 'link':d['link']}
        c = 1
        for row in reader:
            if len(row) >1:
                if row[id_col] in iso_dict:
                    _sample_dict = {'id':iso_dict[row[id_col]]}
                    # if d['link'] == 'snp-distances':
                    #     # print(_sample_dict)
                    for col in row:
                        cols.append(col)
                        _sample_dict[col] = row[col]
                    # if d['link'] == 'snp-distances':
                    #     # print(_sample_dict)
                    tables[d['link']]['table'].append(_sample_dict)
                else:
                    _data_dict = {'id':c}
                    c = c + 1
                    for col in row:
                        cols.append(col)
                        _data_dict[col] = row[col]
                    tables[d['link']]['table'].append(_data_dict)
                    

    cols = list(set(cols))
    if d['columns'] != []:
        columns[d['link']] = d['columns']
    elif 'mlst' in d['link']:
        _c = [{'title':'Isolate','field':'Isolate',"headerFilter":"input","headerFilterPlaceholder":"Search isolate"},{'title':'Scheme','field':'Scheme',"headerFilter":"input","headerFilterPlaceholder":"Search scheme"},{'title':'ST','field':'ST',"headerFilter":"input","headerFilterPlaceholder":"Search ST"}]
        _cls = sorted([c for c in cols if c not in ['Isolate','ST','Scheme']])
        
        if _cls != []:
            for i in _cls:
                _c.append({'title':f"{i}",'field':f"{i}","headerFilter":"input","headerFilterPlaceholder":"Search allele"})
        columns[d['link']] = _c
    
    elif 'snp-distances' in d['link']:
        _c = [{'title':'Isolate','field':'Isolate',"headerFilter":"input","headerFilterPlaceholder":"Search isolate"}]
        _cls = sorted([c for c in cols if c not in ['Isolate']])
        for i in _cls:
            _c.append({'title':f"{i}",'field':f"{i}", "headerVertical":"true", "vertAlign":"bottom"})
        columns[d['link']] = _c
    
    else:
        _c = [{'title':'Isolate','field':'Isolate',"headerFilter":"input","headerFilterPlaceholder":"Search isolate"}]
        _cls = sorted([c for c in cols if c not in ['Isolate']])
        for i in _cls:
            _c.append({'title':f"{i}",'field':f"{i}","headerFilter":"input","headerFilterPlaceholder":f"Search {i}"})
        columns[d['link']] = _c
    
    return tables,columns,comment

def _get_tables(_data, wd, isos):

    # {
    # link: [
    #       {
                # id: some_int, column_1: data, column_2:data
    #        },
    # ]}
    iso_dict = dict(zip(isos,range(1,len(isos)+1))) # for id in table
    
    tables = {}
    columns = {}
    comment = {}
    
    for d in _data:
        # print(d)
        if d['link'] == 'pan-genome':
            id_col = 'Genes'
        elif d['type'] == 'table' or d['type'] == 'matrix':
            id_col = 'Isolate'
        else:
            id_col = ''
        # print(id_col)
        tables,columns,comment = _generate_table(d = d, wd = wd, iso_dict= iso_dict, columns= columns,tables=tables,comment=comment,id_col=id_col)
    
    # print(columns["snp-distances"])
    
    return tables,columns,comment
    # pass

def _compile(args):
    
    # get analysis dict
    _dict = json.load(open(f"{pathlib.Path(args.template_dir, 'report_analysis.json')}", 'r'))
    # print(_dict[args.pipeline])
    # for d in _dict[args.pipeline]:
        # print(d)
    
    isos = _get_isos(wd = args.launchdir, iso_list=args.isolates)
    # print(isos)
    reporthtml = pathlib.Path('report.html')
    # # path to html template
    indexhtml = pathlib.Path(args.template_dir,'index.html') 
    tables,columns,comment = _get_tables(_data = _dict[args.pipeline], wd = args.launchdir, isos = isos)
    version_head,version_body = _get_versions(wd = args.launchdir)
    data = {
        'newick' :'',
        'job_id':args.job_id[0],
        'pipeline':args.pipeline,
        'date':args.day,
        'user':args.user,
        'tables':tables,
        'columns': columns,
        'comment':comment,
        'phylo': 'phylo' if args.iqtree == 'true' else 'no_phylo',
        'num_isos':len(isos),
        'version_head': version_head,
        'version_body':version_body,
        'snp_distances': _plot_distances(wd = args.launchdir),
        'snp_density': _plot_snpdensity(reference=args.reference, wd = args.launchdir, isos = isos),
        'pan_svg': _get_pan_genome(d = _dict[args.pipeline], wd = args.launchdir) if args.pipeline == 'pluspan' else ''
        }
    
    data['newick'] = _get_tree_string(pipeline = args.pipeline, wd = args.launchdir, phylo = args.iqtree)
    
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
        default = '',
        nargs='+')
    parser.add_argument('--isolates',
        help = '',
        default = '')
    parser.add_argument('--reference',
        help = '',
        default = '')
    parser.add_argument('--iqtree',
        help = '',
        default = '')
    
    
    
    
    parser.set_defaults(func=_compile)
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

