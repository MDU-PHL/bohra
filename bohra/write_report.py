# This script will draw a tree from a newick using a recursive approach to extract positions of nodes (adapted for my needs from https://github.com/plotly/dash-phylogeny/blob/dev/app.py) and placement with svgwrite.


import svgwrite, re, sys
from svgwrite import cm, mm
from Bio import Phylo
import jinja2, pathlib, pandas, numpy, re
from bokeh.plotting import figure, show, output_file
from bokeh.embed import components
# from bokeh.io import export_png
# import PyQt5
# from ete3 import Tree, TreeStyle, NodeStyle, TextFace

class Tree:

    def read_treefile(self,filename):
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

    def get_y_coordinates(self,tree, dist=.5):
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
        dwg = svgwrite.Drawing(filename=outpath, debug=True)
        dwg.viewbox(minx=-37.8, miny=(minheight*37.8)-37.8, width=1024, height=((maxheight*37.8)))
        # set horizontal lines
        branches = dwg.add(dwg.g(id='hline', stroke='black'))
        # set vertical lines
        vlines = dwg.add(dwg.g(id='vline', stroke='black'))
        # add tips/nodes
        tips = dwg.add(dwg.g(id='shapes', fill="#3973ac"))
        # add labels
        labels = dwg.add(dwg.g(font_size=8))
        for t in tree_coords:
            if t['type'] == 'horizontal':
                branches.add(dwg.line(start=((t['x0']*f)*cm, t['y0']*cm), end=((t['x1']*f)*cm, t['y1']*cm)))
                
                if t['nodename'] in terms:
                    labels.add(dwg.text(t['nodename'], ((((t['x1']*f) + 0.1)*cm, (t['y1']*cm))),style = 'font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif, "Apple Color Emoji", "Segoe UI Emoji", "Segoe UI Symbol", "Noto Color Emoji"'))
                    circleA = dwg.circle(center=((t['x1']*f)*cm, t['y1']*cm), r='0.05cm')
                    tips.add(circleA)
            elif t['type'] == 'vertical':
                vlines.add(dwg.line(start = ((t['x0']*f)*cm, t['y0']*cm), end = ((t['x1']*f)*cm, t['y1']*cm)))
        dwg.save()



class Report:

    def write_tables(self,reportdir, table):

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
        tablehead = [f"<th>{column}</th>" for column in header]
        # for body seqtablebody
        body = []
        for i in range(1,len(data)):
            row = [f"<tr>"]
            for d in data[i].split('\t'):
                row.append(f"<td>{d}</td>")
            row.append(f"</tr>")
            body = body + row
        return('\n'.join(tablehead),'\n'.join(body))

    def get_table_data(self,reportdir, pipeline, td):
        

        for tabletype in range(len(td)):
            table = td[tabletype]['file']
            td[tabletype]['head'], td[tabletype]['body'] = self.write_tables(reportdir=reportdir, table=table)
        
        return(td)

    def plot_histogram(self,series, xlabel, ylabel, bins, color = "#3973ac", width = 1000, plot_height = 200):
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

        core = reportdir / 'core.tab'
        
        df = pandas.read_csv(core, sep = '\t')
        # get a list of isolate names
        names = list(df.columns[3:len(df.columns)])
        # if the there is no snp in the isolate (ie same as ref then mak na - then easy to drop)
        for i in names:
                df[i]=numpy.where(df['REF'] == df[i], numpy.nan, df[i] )
        # collect positions to get allow for histogram and dropna (no snp)
        melted_df = pandas.melt(df, id_vars=['POS'], value_vars=names)
        melted_df = melted_df.dropna()
        melted_df = melted_df.sort_values(by= ['POS'])
        # generate histogram
        snpdensityscript, snpdensitydiv = self.plot_histogram(series=melted_df['POS']/1000,xlabel="Genome Position (MB)", ylabel="SNPS",bins=10000)
        
        
        # return dictionary
        return(snpdensityscript,snpdensitydiv)

    def plot_distances(self,reportdir):

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
        distancescript, distancediv = self.plot_histogram(series=melted_df['value'], xlabel="SNP distances", ylabel="Frequency", bins=100)
        # td['pairwisedistance'] = {'script':distancecript, 'div':distancediv}
        # return dictionary
        return(distancescript, distancediv)


    def get_tree_image_path(self,reportdir):
        # get tree
        nwk=f"{reportdir / 'core.treefile'}"
        out = f"{reportdir / 'core_tree.svg'}"
        tree = Tree()
        tree.main(treepath=nwk, outpath=out)
        return(f"core_tree.svg")

    def main(self,workdir, resources, job_id, gubbins = False, pipeline = 'sa'):
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
        # table dictionary for each option
        
        td = [{'file':'seqdata.tab', 'title':'Sequence Data', 'link': 'sequencedata', 'type' : 'table'}]
        core_genome_td = {'file': 'core_genome.tab', 'title': 'Core Genome', 'link':'coregenome', 'type':'table'}
        snp_density_td = {'title': 'SNP density', 'link':'snpdensity', 'type':'graph'}
        core_phylogeny_td = {'title': 'Core Phylogeny', 'link':'corephylogeny', 'file': 'core.treefile', 'type': 'tree'}
        snp_distance_td = {'file': 'distances.tab', 'title':'SNP distance matrix', 'type':'matrix', 'link':'distances'}
        # list of snp tasks
        s_td = [core_genome_td,snp_density_td,core_phylogeny_td, snp_distance_td]
        species_id_td = {'file':'species_identification.tab', 'title': 'Species Identification', 'type': 'table', 'link':'speciesid'}
        mlst_td = {'file':'mlst.tab', 'title':'MLST', 'type':'table', 'link':'mlst'}
        resistome_td = {'file':'resistome.tab', 'title':'Resistome', 'type':'table', 'link':'resistome'}
        # list of assembly tasks
        assembly_stat_td = {'file': 'assembly.tab', 'title':'Assembly', 'type':'table', 'link':'assembly'}
        a_td = [assembly_stat_td, species_id_td, mlst_td, resistome_td]
        roary_td = {'file':'summary_statistics.txt', 'title':'Pan Genome', 'type': 'pan', 'image': 'pan_genome.svg', 'link':'pangenome'}
        
        if pipeline == 's':
            td.extend(s_td)
        elif pipeline == 'a':
            td.extend(a_td)
        elif pipeline == 'sa':
            a_td.extend(s_td)
            td.extend(a_td)
            
            # td.extend(s_td)
        else:
            # a_ll = td.extend()
            td.extend(a_td)
            td.extend(s_td)
            td.extend(roary_td)
        # add data to sections
        
        for t in range(len(td)):
            if td[t]['type'] == 'table':
                td[t]['head'], td[t]['body'] = self.write_tables(reportdir=reportdir, table=td[t]['file'])
            if td[t]['type'] == 'tree':
                td[t]['image'] = self.get_tree_image_path(reportdir = reportdir)
            if td[t]['type'] == 'pan':
                td[t]['head'], td[t]['body'] = self.write_tables(reportdir=reportdir, table=td[t]['file'])
            if td[t]['link'] == 'distances':
                td[t]['head'], td[t]['body'] = self.write_tables(reportdir=reportdir, table=td[t]['file'])
                td[t]['script'], td[t]['div'] = self.plot_distances(reportdir=reportdir)
            if td[t]['link'] == 'snpdensity':
                td[t]['script'], td[t]['div'] = self.plot_snpdensity(reportdir= reportdir)
        
        
        # fill template
        report_template = jinja2.Template(pathlib.Path(indexhtml).read_text())
        reporthtml.write_text(report_template.render(td = td, job_id = job_id))
       
        return(True)

if __name__ == '__main__':
    report = Report()
    wd = f"{sys.argv[1]}"
    p = f"{sys.argv[3]}"
    i = f"{sys.argv[4]}"
    report.main(resources=f"{sys.argv[2]}", workdir=wd, pipeline = p, job_id = i)