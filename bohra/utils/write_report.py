
import re, sys, subprocess, toml
import jinja2, pathlib, pandas, numpy, re
from packaging import version
import datetime
from snakemake import shell

class Report:
    '''
    A class to generate the tables and figures for use in the report.html
    '''
    
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
        print(p)
        reporthtml = pathlib.Path('report.html')
        print(reporthtml)
        # path to html template
        indexhtml = pathlib.Path(resources,'index.html') # replace with template
        print(indexhtml)
        print(inputs)
        tml = {}
        with open(inputs, 'r') as t:
            for line in t.readlines():
                if 'td' not in line:
                    d = toml.loads(line)
                    print(len(d))
                    for k in d:
                        print(k)
                        tml[k] = d[k]
            

        # tml = toml.loads(tml_string) 
        # print(tml)
        report_template = jinja2.Template(pathlib.Path(indexhtml).read_text())
        reporthtml.write_text(report_template.render(tml))
    

{input} {params.work_dir} {params.template_path}
inputs = snakemake.input
resources = snakemake.params.template_path
workdir = snakemake.params.work_dir