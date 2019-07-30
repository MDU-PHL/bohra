
import svgwrite, re, sys, subprocess
from svgwrite import cm, mm
from Bio import Phylo
import jinja2, pathlib, pandas, numpy, re
from packaging import version
import datetime
# from bokeh.io import export_png
# import PyQt5
# from ete3 import Tree, TreeStyle, NodeStyle, TextFace

class Tree:
    '''
    This script will draw a tree from a newick using a recursive approach to extract positions of nodes (adapted for my needs from https://github.com/plotly/dash-phylogeny/blob/dev/app.py) and placement with svgwrite.
    '''
    def read_treefile(self,filename):
        '''
        Open the treefile
        '''
        tree = Phylo.read(filename, "newick")
        return tree 

    def get_x_coordinates(self,tree):
        """Associates to  each clade an x-coord.
        returns dict {clade: x-coord}
        """
        xcoords = tree.depths()
        # tree.depth() maps tree clades to depths (by branch length).
        # returns a dict {clade: depth} where clade runs over all Clade instances of the tree, and depth
        # is the distance from root to clade

        #  If there are no branch lengths, assign unit branch lengths
        if not max(xcoords.values()):
            xcoords = tree.depths(unit_branch_lengths=True)
        return xcoords

    def get_y_coordinates(self,tree, dist=.8):
        """
        returns  dict {clade: y-coord}
        The y-coordinates are  (float) multiple of integers (i*dist below)
        dist depends on the number of tree leafs
        """
        maxheight = tree.count_terminals()  # Counts the number of tree leafs.
        # Rows are defined by the tips/leafs
        ycoords = dict((leaf, maxheight - i * dist) for i, leaf in enumerate(reversed(tree.get_terminals())))
        
        def calc_row(clade):
            for subclade in clade:
                if subclade not in ycoords:
                    calc_row(subclade)
            ycoords[clade] = (ycoords[clade.clades[0]] +
                            ycoords[clade.clades[-1]]) / 2

        if tree.root.clades:
            calc_row(tree.root)
        return ycoords

    def get_clade_lines(self,orientation='horizontal',y_curr=0, x_start=0, x_curr=0, y_bot=0, y_top=0):
        """define a shape of type 'line', for branch
        """
        branch_line = dict()
        if orientation == 'horizontal':
            branch_line.update(x0=x_start,
                            y0=y_curr,
                            x1=x_curr,
                            y1=y_curr)
        elif orientation == 'vertical':
            branch_line.update(x0=x_curr,
                            y0=y_bot,
                            x1=x_curr,
                            y1=y_top,
                            type = 'vertical')
        else:
            raise ValueError("Line type can be 'horizontal' or 'vertical'")

        return branch_line


    def draw_clade(self,clade, x_start, line_shapes, line_width=1, x_coords=0, y_coords=0):
        """Recursively draw the tree branches, down from the given clade"""

        x_curr = x_coords[clade]
        y_curr = y_coords[clade]

        # Draw a horizontal line from start to here
        branch_line = self.get_clade_lines(orientation='horizontal', y_curr=y_curr, x_start=x_start, x_curr=x_curr)
        # add in name of node
        branch_line['nodename'] = clade.name
        branch_line['type'] = 'horizontal'
        line_shapes.append(branch_line)
        if clade.clades:
            # Draw a vertical line connecting all children
            y_top = y_coords[clade.clades[0]]
            y_bot = y_coords[clade.clades[-1]]
            line_shapes.append(self.get_clade_lines(orientation='vertical', x_curr=x_curr, y_bot=y_bot, y_top=y_top))

    # Draw descendant
            for child in clade:
                self.draw_clade(child, x_curr, line_shapes, x_coords=x_coords, y_coords=y_coords)
        
        return(line_shapes)

    def main(self,treepath, outpath):
        '''
        Use svgwrite to generate a tree
        '''
        # get the tree using bio phylo
        tree = self.read_treefile(treepath)
        # get the x and y co-ordinates for nodes
        xs = self.get_x_coordinates(tree)
        ys = self.get_y_coordinates(tree)
        terms = [t.name for t in tree.get_terminals()]
        line_shapes = [] # empty list to begin 
        # line shapes is a list of dictionaries for each node and branch
        tree_coords = self.draw_clade(tree.root, 0, line_shapes, x_coords=xs,y_coords=ys)
    
        # find the maximum x position for generation of a figure which will scale
        maxlength = 0
        minheight = len(terms)
        maxheight = 0
        for t in tree_coords:
            if t['type'] == 'horizontal':
                if t['x1'] > float(maxlength):
                    maxlength = t['x1']
                if t['y1'] < float(minheight):
                    minheight = t['y1']
                if t['y1'] > float(maxheight):
                    maxheight = t['y1']
        f = 10/maxlength # factor to multiply the
        # a drawing object
        # viewbox is calculated based on 1cm = 37.8 pixels
        svg_text = [f"<svg baseProfile=\"full\" version=\"1.1\" viewBox=\"-37.8,{(minheight*37.8)-37.8},1000,{(maxheight*37.8)-37.8}\" ><defs />\""]
        
        for t in tree_coords:
            if t['type'] == 'horizontal':
                # branches.add(dwg.line(start=((t['x0']*f)*cm, t['y0']*cm), end=((t['x1']*f)*cm, t['y1']*cm)))
                svg_text.append(f"<line x1=\"{(t['x0']*f)*cm}\" x2=\"{(t['x1']*f)*cm}\" y1=\"{t['y0']*cm}\" y2=\"{t['y1']*cm}\" stroke=\"black\"/>")
                if t['nodename'] in terms:
                    svg_text.append(f"<text class = \"tiplab {t['nodename']}\" x=\"{((t['x1']*f) + 0.1)*cm}\" y=\"{t['y1']*cm}\">{t['nodename']}</text>")
                    # labels.add(dwg.text(t['nodename'], ((((t['x1']*f) + 0.1)*cm, (t['y1']*cm))),style = 'font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif, "Apple Color Emoji", "Segoe UI Emoji", "Segoe UI Symbol", "Noto Color Emoji"'))
                    # circleA = dwg.circle(center=((t['x1']*f)*cm, t['y1']*cm), r='0.05cm')
                    svg_text.append(f"<circle cx=\"{(t['x1']*f)*cm}\" cy=\"{t['y1']*cm}\" r=\"0.05cm\" />")
                    # tips.add(circleA)
                elif (t['nodename'] not in terms) and (t['nodename'] != 'None'):
                    svg_text.append(f"<text class = \"branch-support\" x=\"{((t['x1']*f) + 0.1)*cm}\" y=\"{t['y1']*cm}\" style=\"font-size:smaller; color:#3973ac; display:none;\">{t['nodename']}</text>")
            elif t['type'] == 'vertical':
                # vlines.add(dwg.line(start = ((t['x0']*f)*cm, t['y0']*cm), end = ((t['x1']*f)*cm, t['y1']*cm)))
                svg_text.append(f"<line x1=\"{(t['x0']*f)*cm}\" x2=\"{(t['x1']*f)*cm}\" y1=\"{t['y0']*cm}\" y2=\"{t['y1']*cm}\" stroke = \"black\"/>")
        # dwg.save()
        # svg_text.append('</svg')
        return('\n'.join(svg_text))



class Report:
    '''
    A class to generate the tables and figures for use in the report.html
    '''
    def write_tables(self,reportdir, table):
        '''
        Write a table, given a tab delimited file generate a html string
        '''
        # TODO add class isolate id to <tr>
        # TODO add class distances-isolate to tr if matrix and head-isolate to head td
        path = reportdir / f"{table}"
        data = open(path).readlines()
        # for header
        header = data[0].split('\t')
        if 'mlst' in f"{table}":
            # number of alleles
            length = len(data[1].split('\t'))-2
            header = ['Isolate', 'Scheme', 'ST']
            for l in range(1,length):
                header.append(f"Allele_{l}")
        if 'distances.tab' in table:
            tablehead = [f"<th class='{column}-head'>{column}</th>" for column in header]
        else:
            tablehead = [f"<th>{column}</th>" for column in header]
        # for body seqtablebody
        body = []
        for i in range(1,len(data)):
            raw = data[i].split('\t')
            if 'summary_table.tab' in table:
                row = [f"<tr class='{raw[0]} tiplab'>"]
            elif 'distances.tab' in table:
                row = [f"<tr class='distances-{raw[0]}'>"]
            elif 'assembly.tab' in table:
                row = [f"<tr class='{raw[0]}-assembly'>"]
            elif 'species' in table:
                row = [f"<tr class='{raw[0]}-species-identification'>"]
            elif 'core_genome.tab' in table:
                row = [f"<tr class='{raw[0]}-core-genome'>"]
            elif 'mlst.tab' in table:
                row = [f"<tr class='{raw[0]}-mlst'>"]
            elif 'resistome.tab' in table:
                row = [f"<tr class='{raw[0]}-resistome'>"]
            elif 'seqdata.tab' in table:
                row = [f"<tr class='{raw[0]}-sequence-data'>"]
            else:
                row = [f"<tr>"] # TODO add class isolate id to <tr>
            for d in raw:
                row.append(f"<td align=\"center\">{d}</td>")
            row.append(f"</tr>")
            body = body + row
        return('\n'.join(tablehead),'\n'.join(body))

    def get_table_data(self,reportdir, td):
        '''
        input:
            :reportdir: the directory where report files are kept
            :td: the dictionary for data
        output:
            :td: an updated dictionary
        '''

        for tabletype in range(len(td)):
            table = td[tabletype]['file']
            td[tabletype]['head'], td[tabletype]['body'] = self.write_tables(reportdir=reportdir, table=table)
        
        return(td)

    def plot_histogram(self,series, xlabel, ylabel, bins, color = "#3973ac", width = 1000, plot_height = 200):
        '''
        generic method for generating a histogram from a dataframe using bokeh
        input:
            :series:the column of the dataframe to be used
            :xlabel: label for x-axis
            :ylabel: label for y-axis
            :bins: the size of the histogram bins
            :color: color of the bars
            :width: width of the graph
            :height: height of the graph
        output:
            :script: the javascript string for insert into html
            :div: the html div for inster in html doc
        '''
        # generate histogram
        hist, edges = numpy.histogram(series, density=True, bins=bins)
        # generate bokeh plot
        p = figure(plot_width=width, plot_height=plot_height)
        p.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:], line_color=color, color = color)
        # style
        p.xaxis.axis_label = xlabel
        p.yaxis.axis_label = ylabel
        p.sizing_mode = "scale_width"
        # save as png
        # export_png(p, filename="snpdensity.png")
        # get div and script and and add to dict
        script, div = components(p)
        # show(p)
        return(script,div)

    def plot_snpdensity(self,reportdir):

        '''
        generate a snp-density accross the genome plot - using the core.tab file
        input:
            :reportdir: the directory where files are kept
        out:
            :snpdensityscript: the javascript string for insert into html
            :snpdenstiydiv: the html div for inster in html doc
        '''

        # open fai file
        core = reportdir / 'core.tab'
        
        df = pandas.read_csv(core, sep = '\t')
        # get a list of isolate names
        names = list(df.columns[3:len(df.columns)])
        # if the there is no snp in the isolate (ie same as ref then mak na - then easy to drop)
        for i in names:
                df[i]=numpy.where(df['REF'] == df[i], numpy.nan, df[i])
        # collect positions to get allow for histogram and dropna (no snp)
        melted_df = pandas.melt(df, id_vars=['POS'], value_vars=names)
        melted_df = melted_df.dropna()
        melted_df = melted_df.sort_values(by= ['POS'])
        # generate histogram
        # snpdensityscript, snpdensitydiv = self.plot_histogram(series=melted_df['POS']/1000,xlabel="Genome Position (MB)", ylabel="SNPS",bins=10000)
        
        
        # return dictionary
        return(list(melted_df['POS']/1000))

    def plot_distances(self,reportdir):

        '''
        generate a snp-density plot - using the distacnes.tab file
        input:
            :reportdir: the directory where files are kept
        out:
            :distancescript: the javascript string for insert into html
            :distancesdiv: the html div for inster in html doc
        '''
        distance = reportdir / 'distances.tab'

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


    def get_tree_image(self,reportdir):
        '''
        Generate a tree image from a newick
        input:
            :reportdir: the directory where report files are stored
        output:
            string reporesentation of the path to the tree image
        '''
        # get tree
        nwk=f"{reportdir / 'core.treefile'}"
        out = f"{reportdir / 'core_tree.svg'}"
        tree = Tree()
        
        return(tree.main(treepath=nwk, outpath=out))

    def get_software_versions(self, software):

        '''
        Given the name of the software, find the version
        input:
            :software: the name of the software
        output:
            a string in the form of 'Name_of_Sofware v.X.Y.Z'
        '''

        version_pat = re.compile(r'\bv?(?P<major>[0-9]+)\.(?P<minor>[0-9]+)\.(?P<release>[0-9]+)(?:\.(?P<build>[0-9]+))?\b')

        if software == 'snp-dists':
            v = '-v'
        else:
            v = '--version'
        
        if software in ['snippy', 'prokka']:
            sft = subprocess.run([software, v], stderr=subprocess.PIPE)
            sft = sft.stderr.decode().strip()
        else:
            sft = subprocess.run([software, v], stdout=subprocess.PIPE)
            sft = sft.stdout.decode().strip()
        v = version_pat.search(sft)
        v = v.group()
        sft_version = f"{software} v.{v}"
        return(sft_version)
    
    def make_dict_versions(self, tools):
        '''
        Called by get_software_file to make a dictionary of tools used
        input:
            :tools: a list of tools
        output:
            a dictionary with tools and versions.
        '''
        tool_dict = {}
        for t in tools:
               v = self.get_software_versions(t)
               tool_dict[t] = v
        return(tool_dict)

    def get_software_file(self, reportdir, pipeline, assembler ):
        '''
        get the versions of software on the system at completion of the piepline
        input:
            :reportdir: the directory where report files are stored
            :pipeline: the type of pipeline
            :assembler: the assembler used in the pipeline
        '''
        snippy_tools = ['snippy', 'snippy-core', 'snp-dists', 'iqtree']
        assembly_tools = ['mlst', 'kraken2', 'prokka', 'abricate', assembler]
        
        if pipeline == 's':
            tool_dict = self.make_dict_versions(snippy_tools)
        elif pipeline == 'a':
            tool_dict = self.make_dict_versions(assembly_tools)
        elif pipeline == 'sa':
            snippy_tools.extend(assembly_tools)
            tool_dict = self.make_dict_versions(snippy_tools)
        elif pipeline == 'all':
            snippy_tools.extend(assembly_tools)
            snippy_tools.append('roary')
            tool_dict = self.make_dict_versions(snippy_tools)
        
        
        versions = ['Software versions']
        for t in tool_dict:
            versions.append(tool_dict[t])
        
        p = reportdir / 'software_versions.tab'

        p.write_text('\n'.join(versions))

    def merge_dfs(self,start, added):
        if start.empty:
            start = added
        else:
            start = start.merge(added)
        return(start)

    def generate_summary(self, reportdir):
        '''
        function to generate a summary table
        '''
        tabs = [t for t in reportdir.iterdir() if f"{t.suffix}" == '.tab']
        summary_df = pandas.DataFrame()
        df_list = []
        # print(tabs)
        for tab in tabs:
            # print(tab)
            # print(df)
            if 'species' in f"{tab}":
                species = pandas.read_csv(tab, sep = '\t')
                species = species[['Isolate', '#1 Match']]
                summary_df = self.merge_dfs(summary_df, species)
            elif 'seqdata' in f"{tab}":
                seq = pandas.read_csv(tab, sep = '\t')
                seq = seq[['Isolate', 'Estimated depth']]
                summary_df = self.merge_dfs(summary_df, seq)
            elif 'assembly' in f"{tab}":
                assembly = pandas.read_csv(tab, sep = '\t')
                assembly = assembly[['Isolate', '# Contigs']]
                summary_df = self.merge_dfs(summary_df, assembly)
            elif 'mlst' in f"{tab}":
                mlst = pandas.read_csv(tab, sep = '\t', skiprows=1, header=None)
                mlst = mlst.rename(columns = {0:'Isolate', 2:'ST'})
                mlst = mlst[['Isolate', 'ST']]
                summary_df = self.merge_dfs(summary_df, mlst)
            elif 'core_genome' in f"{tab}":
                core = pandas.read_csv(tab, sep = '\t')
                core = core[['Isolate', '% USED']]
                summary_df = self.merge_dfs(summary_df, core)
        # print(mlst)
        # print(summary_df)
        summary_df = summary_df.rename(columns={'#1 Match': 'Species'})
        summary_file = reportdir / 'summary_table.tab'
        summary_df.to_csv(summary_file, sep = '\t', index = False)
        

    def main(self,workdir, resources, job_id, assembler = 'shovill', gubbins = False, pipeline = 'sa'):
        '''
        main function of the report class ties it all together
        input:
            :workdir: job directory
            :resources: the directory where templates are stored
            :job_id: the id of the job (for html header)
            :assembler: assembler used - for versions
            :gubbins: not in use yet
            :pipeline: the type of pipeline - default is sa = snippy and assembly
        '''

        # set up paths variables
        
        # path to report data
        reportdir = pathlib.Path(workdir,'report')
        reporthtml = reportdir / 'report.html'
        # path to html template
        indexhtml = pathlib.Path(resources,'index.html') # replace with template 
        # copy scc across to report dir
        csstemplate = jinja2.Template(pathlib.Path(resources,'job.css').read_text())
        csstarget = reportdir / 'job.css'
        csstarget.write_text(csstemplate.render())
        # save tool table
        self.get_software_file(reportdir = reportdir, pipeline = pipeline, assembler = assembler)
        
        # table dictionary for each option
        
        td = [{'file':'seqdata.tab', 'title':'Sequence Data', 'link': 'sequence-data', 'type' : 'table'}, {'file':'summary_table.tab','title':'Summary', 'link':'summary', 'type':'summary'}]
        
        # TODO edit links to be title lower case separated by a -
        core_genome_td = {'file': 'core_genome.tab', 'title': 'Core Genome', 'link':'core-genome', 'type':'table'}
        snp_density_td = {'title': 'SNP density', 'link':'snp-density', 'type':'graph'}
        core_phylogeny_td = {'title': 'Phylogeny', 'link':'phylogeny', 'file': 'core.treefile', 'type': 'tree'}
        snp_distance_td = {'file': 'distances.tab', 'title':'SNP distances', 'type':'matrix', 'link':'snp-distances'}
        # list of snp tasks
        s_td = [core_genome_td,snp_density_td,core_phylogeny_td, snp_distance_td]
        species_id_td = {'file':'species_identification.tab', 'title': 'Species Identification', 'type': 'table', 'link':'species-identification'}
        mlst_td = {'file':'mlst.tab', 'title':'MLST', 'type':'table', 'link':'mlst'}
        resistome_td = {'file':'resistome.tab', 'title':'Resistome', 'type':'table', 'link':'resistome'}
        # list of assembly tasks
        assembly_stat_td = {'file': 'assembly.tab', 'title':'Assembly', 'type':'table', 'link':'assembly'}
        a_td = [assembly_stat_td, species_id_td, mlst_td, resistome_td]
        
        # print(roary_td)
        # print(pipeline)
        # print(td)
        if pipeline == 's':
            td.extend(s_td)
            tables =['core-genome', 'snp-distances', 'sequence-data', 'versions']
            modaltables =['core-genome',  'sequence-data']
            display = f"display:inline;"
        elif pipeline == 'a':
            td.extend(a_td)
            tables =['mlst', 'assembly', 'resistome', 'sequence-data','species-identification', 'versions']
            modaltables = tables
            display = f"display:none;"
        elif pipeline == 'sa':
            a_td.extend(s_td)
            td.extend(a_td)
            tables =['core-genome', 'snp-distances', 'mlst', 'assembly', 'resistome', 'sequence-data','species-identification', 'versions'], 
            modaltables = ['core-genome',  'mlst', 'assembly', 'resistome', 'sequence-data', 'species-identification']
            display = f""
            # td.extend(s_td)
        elif pipeline == 'all':
            # a_ll = td.extend()
            roary_td = [{'file':'summary_statistics.txt', 'title':'Pan Genome', 'type': 'pan', 'image': f"{pathlib.Path('pan_genome.svg').open().read()}", 'link':'pan-genome'}]
            td.extend(a_td)
            print(a_td)
            td.extend(s_td)
            print(s_td)
            td.extend(roary_td)
            print(roary_td)
            tables =['core-genome', 'snp-distances', 'mlst', 'assembly', 'resistome', 'sequence-data','species-identification', 'pan-genome', 'versions']
            modaltables = ['core-genome',  'mlst', 'assembly', 'resistome', 'sequence-data', 'species-identification']
            display = f""

        # get versions of software
        versions_td = {'file': 'software_versions.tab', 'title': 'Tools', 'type': 'versions', 'link':'versions'}
        td.append(versions_td)
        snpdistances = ''
        snpdensity = ''
        # add data to sections
        # print(td)
        for t in range(len(td)):
            # print(t)
            # TODO if table add a modal modal + link and link will be title lowercase with hyphen
            if td[t]['type'] == 'table':
                td[t]['head'], td[t]['body'] = self.write_tables(reportdir=reportdir, table=td[t]['file'])
            if td[t]['type'] == 'tree':
                td[t]['image'] = self.get_tree_image(reportdir = reportdir)
            if td[t]['type'] == 'pan':
                td[t]['head'], td[t]['body'] = self.write_tables(reportdir=reportdir, table=td[t]['file'])

            if td[t]['link'] == 'snp-distances':
                td[t]['head'], td[t]['body'] = self.write_tables(reportdir=reportdir, table=td[t]['file'])
                snpdistances = self.plot_distances(reportdir=reportdir)
            if td[t]['link'] == 'snp-density':
                snpdensity= self.plot_snpdensity(reportdir= reportdir)
            if td[t]['type'] == 'versions':
                td[t]['head'], td[t]['body'] = self.write_tables(reportdir=reportdir, table=td[t]['file'])
            if td[t]['type'] == 'summary':
                self.generate_summary(reportdir = reportdir)
                td[t]['head'], td[t]['body'] = self.write_tables(reportdir = reportdir, table = td[t]['file'])
        
        # fill template
        date = datetime.datetime.today().strftime("%d/%m/%y")
        report_template = jinja2.Template(pathlib.Path(indexhtml).read_text())
        reporthtml.write_text(report_template.render(display = display,tables = tables,td = td, job_id = job_id, pipeline = pipeline, snpdistances=snpdistances, snpdensity = snpdensity, modaltables = modaltables, date = date))
    #    TODO pass a list of links for the javascript section called 'table'
    # TODO pass the value of the graphs as separate variable 
        return(True)

if __name__ == '__main__':
    report = Report()
    wd = f"{sys.argv[1]}"
    p = f"{sys.argv[3]}"
    i = f"{sys.argv[4]}"
    a = f"{sys.argv[5]}"
    report.main(resources=f"{sys.argv[2]}", workdir=wd, pipeline = p, job_id = i, assembler=a)