
import re, sys, subprocess, toml
import jinja2, pathlib, pandas, numpy, re
from packaging import version
import datetime

class Report:
    '''
    A class to generate the tables and figures for use in the report.html
    '''
    # def write_tables(self,reportdir, table):
    #     '''
    #     Write a table, given a tab delimited file generate a html string
    #     '''
    #     # TODO add class isolate id to <tr>
    #     # TODO add class distances-isolate to tr if matrix and head-isolate to head td
    #     path = reportdir / f"{table}"
    #     data = open(path).readlines()
    #     # for header
    #     header = data[0].split('\t')
    #     if 'mlst' in f"{table}":
    #         # number of alleles
    #         length = len(data[1].split('\t'))-2
    #         header = ['Isolate', 'Scheme', 'ST']
    #         for l in range(1,length):
    #             header.append(f"Allele_{l}")
    #     if 'distances.tab' in table:
    #         tablehead = [f"<th class='{column}-head'>{column}</th>" for column in header]
    #     else:
    #         tablehead = [f"<th>{column}</th>" for column in header]
    #     # for body seqtablebody
    #     body = []
    #     for i in range(1,len(data)):
    #         raw = data[i].split('\t')
    #         if 'summary_table.tab' in table:
    #             row = [f"<tr class='{raw[0]} tiplab'>"]
    #         elif 'distances.tab' in table:
    #             row = [f"<tr class='distances-{raw[0]}'>"]
    #         elif 'assembly.tab' in table:
    #             row = [f"<tr class='{raw[0]}-assembly'>"]
    #         elif 'species' in table:
    #             row = [f"<tr class='{raw[0]}-species-identification'>"]
    #         elif 'core_genome.tab' in table:
    #             row = [f"<tr class='{raw[0]}-core-genome'>"]
    #         elif 'mlst.tab' in table:
    #             row = [f"<tr class='{raw[0]}-mlst'>"]
    #         elif 'resistome.tab' in table:
    #             row = [f"<tr class='{raw[0]}-resistome'>"]
    #         elif 'seqdata.tab' in table:
    #             row = [f"<tr class='{raw[0]}-sequence-data'>"]
    #         else:
    #             row = [f"<tr>"] # TODO add class isolate id to <tr>
    #         if 'distances.tab' in table:
    #                 for d in range(len(raw)):
    #                     row.append(f"<td align=\"center\" class = \"{raw[0]}_{header[d]}\">{raw[d]}</td>")    
    #             else:
    #                 for d in raw:
    #                     row.append(f"<td align=\"center\">{d}</td>")
    #             row.append(f"</tr>")
    #             body = body + row
    #     return('\n'.join(tablehead),'\n'.join(body))

    

    # def get_table_data(self,reportdir, td):
    #     '''
    #     input:
    #         :reportdir: the directory where report files are kept
    #         :td: the dictionary for data
    #     output:
    #         :td: an updated dictionary
    #     '''

    #     for tabletype in range(len(td)):
    #         table = td[tabletype]['file']
    #         td[tabletype]['head'], td[tabletype]['body'] = self.write_tables(reportdir=reportdir, table=table)
        
    #     return(td)

    # def plot_histogram(self,series, xlabel, ylabel, bins, color = "#3973ac", width = 1000, plot_height = 200):
    #     '''
    #     generic method for generating a histogram from a dataframe using bokeh
    #     input:
    #         :series:the column of the dataframe to be used
    #         :xlabel: label for x-axis
    #         :ylabel: label for y-axis
    #         :bins: the size of the histogram bins
    #         :color: color of the bars
    #         :width: width of the graph
    #         :height: height of the graph
    #     output:
    #         :script: the javascript string for insert into html
    #         :div: the html div for inster in html doc
    #     '''
    #     # generate histogram
    #     hist, edges = numpy.histogram(series, density=True, bins=bins)
    #     # generate bokeh plot
    #     p = figure(plot_width=width, plot_height=plot_height)
    #     p.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:], line_color=color, color = color)
    #     # style
    #     p.xaxis.axis_label = xlabel
    #     p.yaxis.axis_label = ylabel
    #     p.sizing_mode = "scale_width"
    #     # save as png
    #     # export_png(p, filename="snpdensity.png")
    #     # get div and script and and add to dict
    #     script, div = components(p)
    #     # show(p)
    #     return(script,div)

    # def plot_snpdensity(self,reportdir, workdir):

    #     '''
    #     generate a snp-density accross the genome plot - using the core.tab file
    #     input:
    #         :reportdir: the directory where files are kept
    #     out:
    #         :snpdensityscript: the javascript string for insert into html
    #         :snpdenstiydiv: the html div for inster in html doc
    #     '''
    #     # helper functions for getting the data into the right format.
    #     def adjust_offset(row, d):
    #         if row['CHR'] in d:
    #             return(int(d[row['CHR']]) + int(row['POS']))
    #         else:
    #             return(int(row['POS']))

    #     def generate_dict(idx_file):
    #         d = {}
    #         pos = 0
    #         offset = 0
    #         with open(f"{idx_file}") as f:
    #             for line in f:
    #                 l = line.split()
    #                 d[l[0]] = offset
    #                 offset = offset + int(l[1])
    #         return(d)
        
    #     # open fai file and generate the dictionary
    #     idx = pathlib.Path(workdir ,'ref.fa.fai')
    #     d = generate_dict(idx)
    #     # get the core file
    #     core = reportdir / 'core.tab'
    #     df = pandas.read_csv(core, sep = '\t')
    #     # get a list of isolate names
    #     names = list(df.columns[3:len(df.columns)])
    #     # if the there is no snp in the isolate (ie same as ref then mak na - then easy to drop)
    #     for i in names:
    #             df[i]=numpy.where(df['REF'] == df[i], numpy.nan, df[i])
    #     # generate the offset value
    #     df['POS_OFFSET'] = df[['CHR', 'POS']].apply(lambda x:adjust_offset(x,d), axis = 1)
    #     # collect positions to get allow for histogram and dropna (no snp)
    #     melted_df = pandas.melt(df, id_vars=['POS_OFFSET'], value_vars=names)
    #     melted_df = melted_df.dropna()
    #     melted_df = melted_df.sort_values(by= ['POS_OFFSET'])
    #     # generate histogram
    #     # snpdensityscript, snpdensitydiv = self.plot_histogram(series=melted_df['POS']/1000,xlabel="Genome Position (MB)", ylabel="SNPS",bins=10000)
    #     contig_breaks = [d[value] for value in d]
        
    #     # return dictionary
    #     return(list(melted_df['POS_OFFSET']/1000))

    # def plot_distances(self,reportdir):

    #     '''
    #     generate a snp-density plot - using the distacnes.tab file
    #     input:
    #         :reportdir: the directory where files are kept
    #     out:
    #         :distancescript: the javascript string for insert into html
    #         :distancesdiv: the html div for inster in html doc
    #     '''
    #     distance = reportdir / 'distances.tab'

    #     df = pandas.read_csv(distance, sep = '\t')
    #     # get a list of isolate names
    #     names = list(df.columns[1:len(df.columns)])
    #     col1 = df.columns[0]
        
    #     # if the there is no snp in the isolate (ie same as ref then mak na - then easy to drop)
        
    #     # collect positions to get allow for histogram and dropna (no snp)
    #     melted_df = pandas.melt(df, id_vars=[col1], value_vars=names)
    #     melted_df = melted_df[melted_df[col1]!= melted_df['variable']]
    #     # generate histogram
    #     # distancescript, distancediv = self.plot_histogram(series=melted_df['value'], xlabel="SNP distances", ylabel="Frequency", bins=100)
    #     # td['pairwisedistance'] = {'script':distancecript, 'div':distancediv}
    #     # return dictionary
    #     return(list(melted_df['value']))


    # def get_tree_string(self,reportdir, tree_file):
    #     '''
    #     Generate a tree image from a newick
    #     input:
    #         :reportdir: the directory where report files are stored
    #     output:
    #         string reporesentation of the path to the tree image
    #     '''
    #     with open(f"{reportdir / tree_file}", 'r') as t:
    #         tree = t.read().strip()

    #     return tree

    # def get_software_versions(self, software):

    #     '''
    #     Given the name of the software, find the version
    #     input:
    #         :software: the name of the software
    #     output:
    #         a string in the form of 'Name_of_Sofware v.X.Y.Z'
    #     '''

    #     version_pat = re.compile(r'\bv?(?P<major>[0-9]+)\.(?P<minor>[0-9]+)\.(?P<release>[0-9]+)(?:\.(?P<build>[0-9]+))?\b')

    #     if software == 'snp-dists':
    #         v = '-v'
    #     else:
    #         v = '--version'
        
    #     if software in ['snippy', 'prokka']:
    #         sft = subprocess.run([software, v], stderr=subprocess.PIPE)
    #         sft = sft.stderr.decode().strip()
    #     else:
    #         sft = subprocess.run([software, v], stdout=subprocess.PIPE)
    #         sft = sft.stdout.decode().strip()
    #     v = version_pat.search(sft)
    #     v = v.group()
    #     sft_version = f"{software} v.{v}"
    #     return(sft_version)
    
    # def make_dict_versions(self, tools):
    #     '''
    #     Called by get_software_file to make a dictionary of tools used
    #     input:
    #         :tools: a list of tools
    #     output:
    #         a dictionary with tools and versions.
    #     '''
    #     tool_dict = {}
    #     for t in tools:
    #            v = self.get_software_versions(t)
    #            tool_dict[t] = v
    #     return(tool_dict)

    # def get_software_file(self, reportdir, pipeline, assembler ):
    #     '''
    #     get the versions of software on the system at completion of the piepline
    #     input:
    #         :reportdir: the directory where report files are stored
    #         :pipeline: the type of pipeline
    #         :assembler: the assembler used in the pipeline
    #     '''
    #     snippy_tools = ['snippy', 'snippy-core', 'snp-dists', 'iqtree']
    #     assembly_tools = ['mlst', 'kraken2', 'prokka', 'abricate', assembler]
        
    #     if pipeline == 's':
    #         tool_dict = self.make_dict_versions(snippy_tools)
    #     elif pipeline == 'a':
    #         tool_dict = self.make_dict_versions(assembly_tools)
    #     elif pipeline == 'sa':
    #         snippy_tools.extend(assembly_tools)
    #         tool_dict = self.make_dict_versions(snippy_tools)
    #     elif pipeline == 'all':
    #         snippy_tools.extend(assembly_tools)
    #         snippy_tools.append('roary')
    #         tool_dict = self.make_dict_versions(snippy_tools)
        
        
    #     versions = ['Software versions']
    #     for t in tool_dict:
    #         versions.append(tool_dict[t])
        
    #     p = reportdir / 'software_versions.tab'

    #     p.write_text('\n'.join(versions))

    # def merge_dfs(self,start, added):
    #     if start.empty:
    #         start = added
    #     else:
    #         start = start.merge(added)
    #     return(start)

    # def generate_summary(self, reportdir):
    #     '''
    #     function to generate a summary table
    #     '''
    #     tabs = [t for t in reportdir.iterdir() if f"{t.suffix}" == '.tab']
    #     summary_df = pandas.DataFrame()
    #     df_list = []
    #     # print(tabs)
    #     for tab in tabs:
    #         # print(df)
    #         if 'species' in f"{tab.name}":
    #             species = pandas.read_csv(tab, sep = '\t')
    #             species = species[['Isolate', '#1 Match']]
    #             summary_df = self.merge_dfs(summary_df, species)
    #             summary_df = summary_df.rename(columns={'#1 Match': 'Species'})
    #         elif 'seqdata' in f"{tab.name}":
    #             seq = pandas.read_csv(tab, sep = '\t')
    #             seq = seq[['Isolate', 'Estimated depth']]
    #             summary_df = self.merge_dfs(summary_df, seq)
    #         elif 'assembly' in f"{tab.name}":
    #             assembly = pandas.read_csv(tab, sep = '\t')
    #             assembly = assembly[['Isolate', '# Contigs']]
    #             summary_df = self.merge_dfs(summary_df, assembly)
    #         elif 'mlst' in f"{tab.name}":
    #             mlst = pandas.read_csv(tab, sep = '\t', skiprows=1, header=None)
    #             mlst = mlst.rename(columns = {0:'Isolate', 2:'ST'})
    #             mlst = mlst[['Isolate', 'ST']]
    #             summary_df = self.merge_dfs(summary_df, mlst)
    #         elif 'core_genome' in f"{tab.name}":
    #             core = pandas.read_csv(tab, sep = '\t')
    #             core = core[['Isolate', '% USED']]
    #             summary_df = self.merge_dfs(summary_df, core)
    #     # print(mlst)
    #     # print(summary_df)
        
    #     summary_file = 'summary_table.tab'
    #     summary_df.to_csv(summary_file, sep = '\t', index = False)
    
    # def get_all(self, pipeline):
    #     '''
    #     return a dictionary for report
    #     '''

    #     td = [{'file':'seqdata.tab', 'title':'Sequence Data', 'link': 'sequence-data', 'type' : 'table'}, {'file':'species_identification.tab', 'title': 'Species Identification', 'type': 'table', 'link':'species-identification'}{'file':'summary_table.tab','title':'Summary', 'link':'summary', 'type':'summary'}]
        
    #     # TODO edit links to be title lower case separated by a -
    #     core_genome_td = {'file': 'core_genome.tab', 'title': 'Core Genome', 'link':'core-genome', 'type':'table'}
    #     snp_density_td = {'title': 'SNP density', 'link':'snp-density', 'type':'graph'}
    #     core_phylogeny_td = {'title': 'Phylogeny', 'link':'phylogeny', 'file': 'core.treefile', 'type': 'tree'}
    #     snp_distance_td = {'file': 'distances.tab', 'title':'SNP distances', 'type':'matrix', 'link':'snp-distances'}
    #     # list of snp tasks
    #     s_td = [core_genome_td,snp_density_td,core_phylogeny_td, snp_distance_td]
    #     mlst_td = {'file':'mlst.tab', 'title':'MLST', 'type':'table', 'link':'mlst'}
    #     resistome_td = {'file':'resistome.tab', 'title':'Resistome', 'type':'table', 'link':'resistome'}
    #     # list of assembly tasks
    #     assembly_stat_td = {'file': 'assembly.tab', 'title':'Assembly', 'type':'table', 'link':'assembly'}
    #     a_td = [assembly_stat_td, species_id_td, mlst_td, resistome_td]
        
        
    #     if pipeline == 's':
    #         td.extend(s_td)
    #         tables =['core-genome', 'snp-distances', 'sequence-data']
    #         modaltables =['core-genome',  'sequence-data','species-identification']
    #         display = f"display:inline;"
    #     elif pipeline == 'a':
    #         td.extend(a_td)
    #         tables =['mlst', 'assembly', 'resistome', 'sequence-data']
    #         modaltables = ['mlst', 'assembly', 'resistome', 'sequence-data','species-identification']
    #         display = f"display:none;"
    #     elif pipeline == 'sa':
    #         a_td.extend(s_td)
    #         td.extend(a_td)
    #         tables =['core-genome', 'snp-distances', 'mlst', 'assembly', 'resistome', 'sequence-data']
    #         modaltables = ['core-genome',  'mlst', 'assembly', 'resistome', 'sequence-data','species-identification']
    #         display = f""
    #         # td.extend(s_td)
    #     elif pipeline == 'all':
    #         # a_ll = td.extend()
            
    #         td.extend(a_td)
    #         td.extend(s_td)
    #         td.extend(roary_td)
    #         tables =['core-genome', 'species-identification','snp-distances', 'mlst', 'assembly', 'resistome', 'sequence-data', 'pan-genome']
    #         modaltables = ['core-genome',  'mlst', 'assembly', 'resistome', 'sequence-data']
    #         display = f""
        
    #     tables.append('versions')
    #     # get versions of software
    #     versions_td = {'file': 'software_versions.tab', 'title': 'Tools', 'type': 'versions', 'link':'versions'}
    #     td.append(versions_td)
    #     snpdistances = ''
    #     snpdensity = ''

    #     return td, snpdistances, snpdensity

    # def get_preview(self):

    #     td = [{'title': 'Preview', 'link':'preview-tree', 'file': 'preview.newick', 'type': 'tree'},
    #            {'file': 'preview_distances.tab', 'title':'Mash distances', 'type':'matrix', 'link':'mash-distances'}
    #     ]
    

    def open_toml(self,tml):

        data = toml.load(tml)

        return data
        
    def write_toml(self, data, output):
        
        with open(output, 'wt') as f:
            toml.dump(data, f)

    def main(self, inputs, workdir, resources):
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
        p = pathlib.Path('.')
        # print(sorted(p.glob('*.fa*')))
        # path to report data
        # reportdir = pathlib.Path(workdir)
        reporthtml = pathlib.Path('report.html')
        # path to html template
        indexhtml = pathlib.Path(resources,'index.html') # replace with template
        # print(indexhtml)
        tml =  self.open_toml(inputs)
        # print(type(tml))
        # print(tml)
        # copy scc across to report dir
        # csstemplate = jinja2.Template(pathlib.Path(resources,'job.css').read_text())
        # csstarget = reportdir / 'job.css'
        # csstarget.write_text(csstemplate.render())
        # newick string
        # tree_file = 'preview.newick' if preview else 'core.treefile'
        # newick = self.get_tree_string(reportdir = reportdir, tree_file = tree_file)
        # # save tool table
        # self.get_software_file(reportdir = reportdir, pipeline = pipeline, assembler = assembler)
        
        # # table dictionary for each option
        # if preview:
        #     td = self.get_preview()
        # else:
        #     td = self.get_all(pipeline = pipeline)
        
        # # add data to sections
        # # print(td)
        # for t in range(len(td)):
        #     # print(t)
        #     # TODO if table add a modal modal + link and link will be title lowercase with hyphen
        #     if td[t]['type'] == 'table':
        #         td[t]['head'], td[t]['body'] = self.write_tables(reportdir=reportdir, table=td[t]['file'])
        #     if td[t]['type'] == 'tree':
        #         td[t]['image'] = self.get_tree_image(reportdir = reportdir)
        #     if td[t]['type'] == 'pan':
        #         td[t]['head'], td[t]['body'] = self.write_tables(reportdir=reportdir, table=td[t]['file'])

        #     if td[t]['link'] == 'snp-distances':
        #         td[t]['head'], td[t]['body'] = self.write_tables(reportdir=reportdir, table=td[t]['file'])
        #         snpdistances = self.plot_distances(reportdir=reportdir)
        #     if td[t]['link'] == 'snp-density':
        #         snpdensity= self.plot_snpdensity(reportdir= reportdir, workdir=workdir)
        #     if td[t]['type'] == 'versions':
        #         td[t]['head'], td[t]['body'] = self.write_tables(reportdir=reportdir, table=td[t]['file'])
        #     if td[t]['type'] == 'summary':
        #         self.generate_summary(reportdir = reportdir)
        #         td[t]['head'], td[t]['body'] = self.write_tables(reportdir = reportdir, table = td[t]['file'])
        
        # # fill template
        # {'file': 'software_versions.tab', 'title': 'Tools', 'type': 'versions', 'link':'versions'}
        report_template = jinja2.Template(pathlib.Path(indexhtml).read_text())
        # print(report_template)
        # print(report_template.render(tml))
        reporthtml.write_text(report_template.render(tml))
    #    TODO pass a list of links for the javascript section called 'table'
    # TODO pass the value of the graphs as separate variable 
        return(True)

if __name__ == '__main__':
    report = Report()
    
    report.main(inputs = f"{sys.argv[1]}", resources=f"{sys.argv[3]}", workdir=f"{sys.argv[2]}")