#!/usr/bin/env python3
import pathlib, subprocess, argparse, datetime, pandas, re, numpy, jinja2, json, csv
import altair as alt
from Bio import SeqIO

alt.data_transformers.disable_max_rows()

def _get_tree_string(pipeline, wd,phylo, report_outdir='report'):
    '''
    Generate a tree image from a newick
    input:
        :reportdir: the directory where report files are stored
    output:
        string reporesentation of the path to the tree image
    '''
    tree_file = pathlib.Path(wd, report_outdir,'preview.newick') if pipeline == 'preview' else pathlib.Path(wd, report_outdir,'core.newick')
    if tree_file.exists() and phylo == 'true':
        with open(f"{tree_file}", 'r') as t:
            tree = t.read().strip()
    else:
        tree = 'No tree available'
    # print(tree)
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

def get_bin_size(_dict):

    sum_len = 0
    for d in _dict:
        sum_len = sum_len + _dict[d]['length']
    
    _maxbins = int(sum_len/3000)
    if _maxbins == 0:
        print(f"Something has gone wrong - the maxbins value should be > 0.")
    return _maxbins

def check_masked(mask_file, df,wd, _dict ):

    masked = []
    if mask_file != '' and pathlib.Path(wd,mask_file).exists():
        print('blocking out masked regions')
        mask = pandas.read_csv(f"{pathlib.Path(wd,mask_file)}", sep = '\t', header = None, names = ['CHR','Pos1','Pos2'])
        mask['CHR'] = mask['CHR'].astype(str)
        for row in mask.iterrows():
            # print(row[1])
            off = _dict[row[1]['CHR']]['offset']
            l = list(range(row[1]['Pos1'] + off, row[1]['Pos2']+off +1))
            masked.extend(l)

    df['mask'] = numpy.where(df['index'].isin(masked), 'masked', 'unmasked')
    
    return df

def get_contig_breaks(_dict):
    for_contigs = []
    for chromosome in _dict:
        if _dict[chromosome]['length'] > 5000:
            for_contigs.append(_dict[chromosome]['length'] + _dict[chromosome]['offset'])

    return for_contigs


def _plot_snpdensity(reference,wd, isos, mask_file = ''):

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
    chromosomes = list(_dict.keys())
    
    _maxbins = get_bin_size(_dict = _dict)
    # collate all snps in snps.tab
    vars = {}
    for i in isos:
        # open snps.tab
        snps = pathlib.Path(wd, i, 'snps.tab')
        if snps.exists():
            tab_file = pandas.read_csv(f"{snps}", dtype = str,sep = '\t')
            for chromosome  in chromosomes:
            # print(var)
            # vars = []
                if chromosome not in vars:
                    vars[chromosome] = {}
                chr = tab_file[tab_file['CHROM'] == chromosome]
                # print(chr)
                if not chr.empty:
                    for row in chr.iterrows():
                        # print(row[1]['POS'])
                        pos = int(row[1]['POS'])
                        if pos not in vars[chromosome]:
                            vars[chromosome][pos] = 1
                        else:
                            vars[chromosome][pos] = vars[chromosome][pos] + 1
    # now generate list for x and y value in graph
    data = {}
    for var in vars:
        for pos in vars[var]:
            offset = _dict[var]['offset']
            data[pos + offset] = vars[var][pos]
    df = pandas.DataFrame.from_dict(data, orient='index',columns=['vars']).reset_index()
    # check if mask file used - if yes grey it out in the graph.
    df = check_masked(mask_file = mask_file, df = df, wd = wd, _dict = _dict)
    # get positions of the contig breaks
    for_contigs = get_contig_breaks(_dict = _dict)
    # set colours
    domain = ['masked', 'unmasked']
    range_ = ['#d9dcde', '#216cb8']
    # do bar graphs
    print(mask_file)
    print(df[df['mask'] == 'masked'])
    # if mask_file != 'no_mask':
    bar = alt.Chart(df).mark_bar().encode(
        x=alt.X('index:Q', bin=alt.Bin(maxbins=_maxbins), title = "Core genome position.", axis=alt.Axis(ticks=False)),
        y=alt.Y('sum(vars):Q',title = "Variants observed (per 500 bp)"),
        color=alt.Color('mask', scale = alt.Scale(domain=domain, range=range_), legend=None)
    )

    # generate list of graphs for addition of vertical lines
    graphs = [bar]
    if for_contigs != []:
        for line in for_contigs:
            graphs.append(alt.Chart().mark_rule(strokeDash=[3, 3], size=1, color = 'grey').encode(x = alt.datum(line)))
        
    chart = alt.layer(*graphs).configure_axis(
                    grid=False
                    ).properties(
                        width = 1200
                    ).interactive()

    chart = chart.to_json()
    
    return chart

def _plot_distances(wd,report_outdir='report'):

    '''
    generate a snp-density plot - using the distacnes.tab file
    input:
        :reportdir: the directory where files are kept
    out:
        :distancescript: the javascript string for insert into html
        :distancesdiv: the html div for inster in html doc
    '''
    distance = f"{pathlib.Path(wd, report_outdir,'distances.tab')}"
    try:
        df = pandas.read_csv(distance, sep = '\t')
        # get a list of isolate names
        names = list(df.columns[1:len(df.columns)])
        col1 = df.columns[0]
        
        # if the there is no snp in the isolate (ie same as ref then mak na - then easy to drop)
        
        # collect positions to get allow for histogram and dropna (no snp)
        melted_df = pandas.melt(df, id_vars=[col1], value_vars=names)
        melted_df = melted_df[melted_df[col1]!= melted_df['variable']]
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

def _subset_versions( _typers, wd, _types, report_outdir='report'):

    if _types != set():
        vs = pandas.read_csv(f"{pathlib.Path(wd, report_outdir, 'versions.txt')}", sep = '\t')
        print(vs)
        to_remove = [i for i in vs['tool'].unique() if (i in _typers['typers'] and (i not in _types))]
        vs = vs[~vs['tool'].isin(to_remove)]
        print(vs)
        vs.to_csv(f"{pathlib.Path(wd, report_outdir, 'versions.txt')}", sep = '\t', index = False)

def _extract_typer( _typers, wd ):
    
    typers = sorted(pathlib.Path(wd).glob(f"*/typer*.txt"))
    print(typers)
    _types = set()
    if typers != []:
        for typer in typers:
            t = typer.name.split('_')[-1].split('.')[0]
            _types.add(t)
    
    print(_types)
    _subset_versions(_typers = _typers, wd = wd, _types = _types)

    return _types


def _get_pan_genome(d, wd, report_outdir='report'):

    image = d[1]['image']
    path = pathlib.Path(wd, report_outdir, image)
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

def _generate_table(d, columns,comment, tables, wd, iso_dict, id_col, report_outdir='report'):
    
    if id_col == '':
        return tables,columns,comment
    cols = []
    try:
        with open(f"{pathlib.Path(wd, report_outdir, d['file'])}", 'r') as f:
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
    except FileNotFoundError:
        print(f"No file found for {d['title']}")

                        

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
        elif d['link'] == 'software-versions':
            id_col = 'tool'
        elif d['link'] == 'pipeline-details':
            id_col = 'detail'
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
    _typers = json.load(open(f"{pathlib.Path(args.template_dir, 'typers.json')}", "r"))
    print(_typers)
    _extract_typer(_typers = _typers, wd = args.launchdir)
    # print(_dict[args.pipeline])
    # for d in _dict[args.pipeline]:
        # print(d)
    
    isos = _get_isos(wd = args.launchdir, iso_list=args.isolates)
    # print(isos)
    reporthtml = pathlib.Path(f'report_{args.pipeline}.html')
    # # path to html template
    indexhtml = pathlib.Path(args.template_dir,'index.html') 
    tables,columns,comment = _get_tables(_data = _dict[args.pipeline], wd = args.launchdir, isos = isos)
    
    data = {
        'newick' :'',
        'job_id':args.job_id[0],
        'pipeline':args.pipeline,
        'date':args.day,
        'user':args.user,
        # 'reference_string':_get_reference_string(reference = args.reference, mask = args.mask, pipeline = args.pipeline),
        'tables':tables,
        'columns': columns,
        'comment':comment,
        'phylo': 'phylo' if args.iqtree == 'true' else 'no_phylo',
        'num_isos':len(isos),
        # 'version_head': version_head,
        # 'version_body':version_body,
        'snp_distances': _plot_distances(wd = args.launchdir, report_outdir=args.report_outdir) if args.pipeline not in ['preview', 'assemble','amr_typing'] else {0:0},
        'snp_density': _plot_snpdensity(reference=args.reference, wd = args.launchdir, isos = isos, mask_file = args.mask, report_outdir=args.report_outdir) if args.pipeline not in ['preview', 'assemble','amr_typing'] else {0:0},
        'pan_svg': _get_pan_genome(d = _dict[args.pipeline], wd = args.launchdir, report_outdir=args.report_outdir) if args.pipeline == 'full' else ''
        }
    
    data['newick'] = _get_tree_string(pipeline = args.pipeline, wd = args.launchdir, phylo = args.iqtree, report_outdir=args.report_outdir)
    
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
    parser.add_argument('--contigs',
        help='',
        default = '')
    parser.add_argument('--mask',
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
    parser.add_argument('--report_outdir',
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

