import pathlib
import os, getpass, shutil, re, psutil
import pandas
import jinja2
import sh
import logging
import filecmp
import datetime
import numpy
import itertools
import subprocess
from Bio import SeqIO, Phylo
from packaging import version
from bohra.utils.write_snakemake import MakeWorkflow
# from bohra.utils.write_report import Report


class RunSnpDetection(object):
    '''
    A class for Bohra pipeline object
    '''
    
    def __init__(self, args):
        # get date and time
        self.now = datetime.datetime.today().strftime("%d_%m_%y_%H")
        self.day = datetime.datetime.today().strftime("%d_%m_%y")
        # get the working directory
        self.workdir = pathlib.Path(args.workdir)
        # path to pipeline resources
        self.resources = pathlib.Path(args.resources)
        # path to reference and mask
        self.ref = pathlib.Path(args.reference)
        # 
        # (args.mask)
        if args.mask:
            rerun_core = self.check_mask(args.mask)
        else:
            self.mask = ''
        # path to input file
        if args.input_file == '':
            self.log_messages('warning', 'Input file can not be empty, please set -i path_to_input to try again')
            raise SystemExit()
        else:
            self.input_file = pathlib.Path(args.input_file)
        # # path to a source file for tracking of reference, jobid and mask
        # self.source_log_path = pathlib.Path(self.workdir, 'source.log')
        # job id
        self.job_id = self._name_exists(args.job_id)
        # other variables
        # min aln 
        self.minaln = args.minaln
        # set threads
        
        # user
        self.user = getpass.getuser()
        # gubbins TODO add back in later!!
        # if not args.gubbins:
        #     self.gubbins = numpy.nan
        # else:
        #     self.gubbins = args.gubbins
        
        self.gubbins = numpy.nan

        if isinstance(args.prefillpath, str):
            self.prefillpath = args.prefillpath
        elif args.mdu:
            self.prefillpath = '/home/seq/MDU/QC'
        else:
            self.prefillpath = ''
        self.force = args.force
        self.dryrun = args.dryrun
        self.pipeline = args.pipeline
        self.cpus = args.cpus
        
        self.assembler = args.assembler
        self.snippy_version = ''
        self.assembler_dict = {'shovill': 'shovill', 'skesa':'skesa','spades':'spades.py'}
        self.set_snakemake_jobs()

    def set_snakemake_jobs(self):
        '''
        set the number of jobs to run in parallel based on the number of cpus from args
        '''
        if int(self.cpus) < int(psutil.cpu_count()):
            self.jobs =  self.cpus
        else:
            self.jobs = 1

    def log_messages(self, type, message):
        '''
        Will log messages to the screen and also add them to job.log
        input:
            :type: type of message info or warning
            :message: message to log

        '''

        if type == 'warning':
            logging.warning(message)
            print(f"WARNING: {message}")
        if type == 'info':
            logging.info(message)
            # print(f"{message}")
    
    def force_overwrite(self):
        '''
        will force pipeline to run in an existing folder - removes isolate and source logs
        '''
        isolatelog = self.workdir / f"isolates.log"
        sourcelog = self.workdir / f"source.log"
        # joblog = self.workdir / f"job.log"
        if isolatelog.exists():
            isolatelog.unlink()
        if sourcelog.exists():
            sourcelog.unlink()
        # if joblog.exists():
        #     joblog.unlink()

        return False

    def check_setup_files(self):
        '''
        check that the working directory, resources directory and the input file exist
        '''
        self.path_exists(self.workdir, v = False) 
        self.path_exists(self.resources, v = False)
        self.path_exists(self.input_file)
    
    def check_snippy(self):
        '''
        check for snippy
        '''
        version_pat = re.compile(r'\bv?(?P<major>[0-9]+)\.(?P<minor>[0-9]+)\.(?P<release>[0-9]+)(?:\.(?P<build>[0-9]+))?\b')
        try:
            snippy = subprocess.run(['snippy', '--version'], stderr=subprocess.PIPE)
            snippy = snippy.stderr.decode().strip()
            self.snippy_version = version_pat.search(snippy)
            self.log_messages('info', f"Snippy {snippy} found. Good job!")
            
            return(version_pat.search(snippy))
        except FileNotFoundError:
            self.log_messages('warning', f"snippy is not installed.")
            raise SystemExit
    
    def check_snippycore(self):
        '''
        check for snippy-core
        '''
        self.check_installation('snippy-core')

    def check_snpdists(self):
        '''
        check for snp-dists
        '''
        self.check_installation('snp-dists')



    def check_iqtree(self):
        '''
        check iqtree
        '''
        self.check_installation('iqtree')

    def check_installation(self,software):
        '''
        Check that software is installed
        input:
            :software: the name of the software - must be command line name 
        '''

        if shutil.which(software):
            self.log_messages('info', f"{software} is installed")
        else:
            self.log_messages('warning', f"{software} is not installed, please check dependencies and try again.")
            raise SystemExit


    def check_assembler(self):
        '''
        check version of assembler
        '''
        ret = 0
        
        self.check_installation(self.assembler_dict[self.assembler])


    def check_assemble_accesories(self):
        '''
        check the assembly accessories mlst, kraken, abricate and prokka
        '''
        accessories = ['mlst', 'kraken2', 'abricate', 'prokka']
        
        for a in accessories:
            self.check_installation(a)
    
    def check_roary(self):
        '''
        check roary is installed
        '''

        self.check_installation('roary')


    def check_deps(self):
        '''
        check dependencies Snippy, snippy-core, snp-dists, iqtree
        '''
        # TODO check all software tools used and is there a way to check database last update??
        # TODO check assemblers
        self.log_messages('info', f"Checking software dependencies")
        if self.pipeline != 'a':
            self.check_snippycore()
            self.check_snpdists()
            self.check_iqtree()
            return(self.check_snippy())
        if self.pipeline != 's':
            self.check_assembler()
            self.check_assemble_accesories()
        if self.pipeline == 'all':
            self.check_roary()
        

    # def check_validation(self, validation_type):
    #     '''
    #     Check if this is a validation run 
    #     '''
    #         if validation_type == 'both':
    #             validate = ['snps', 'clusters']
    #         else:
    #             if validation_type in ['snps', 'clusters']:
    #                 validate = validation_type.split(',')
    #             else:
    #                 self.log_messages('warning', f"{validation_type} is not a valid input. Correct options can be taken from {' '.join(['snps', 'clusters', 'both'])}. Please try again")
    #                 raise SystemExit

    #         return(validate)
    
    

    def run_checks(self):
        '''
        Run checks prior to start - checking all software is installed, if this is a rerun and the input files
        '''
        self.check_setup_files()
        self.check_rerun()
        self.check_deps()
        # check reference
        if self.pipeline != 'a':
            if self.ref == '':
                self.log_messages('warning', f"You are trying call SNPs, a reference file is required. Please try again using '-r path to reference'")
                raise SystemExit
            else:
                self.ref = self.link_file(self.ref)
        
            

    def set_source_log(self):
        '''
        set the reference, mask and id for tracking and potential.
        
            
        '''   
        new_df = pandas.DataFrame({'JobID':self.job_id, 'Reference':f"{self.ref}",'Mask':f"{self.mask}", 
                                    'MinAln':self.minaln, 'Pipeline': self.pipeline, 'CPUS': self.cpus, 'Assembler':self.assembler,
                                    'Gubbins': self.gubbins, 'Date':self.day, 'User':self.user, 'snippy_version':self.snippy_version, 'input_file':f"{self.input_file}",'prefillpath': self.prefillpath}, 
                                    index=[0], )

        source_path = self.workdir / 'source.log'
        if source_path.exists():
            source_df = pandas.read_csv(source_path, '\t')
            source_df = source_df.append(new_df)
        else:
            source_df = new_df
        
        source_df.to_csv(source_path , index=False, sep = '\t')

    def check_rerun(self):
        '''
        CHeck if the job is a rerun of an existing job, if so print message informing user and exit

        '''

        source_path = self.workdir / 'source.log'
        # if the path is a string convert to Path
        if isinstance(source_path, str):
            source_path = pathlib.Path(source_path)
        if source_path.exists():
            self.log_messages('warning',f"This may be a re-run of an existing job. Please try again using rerun instead of detect OR use -f to force an overwrite of the existing job.")
            self.log_messages('warning',f"Exiting....")
            self.log_messages('info',f"{60 * '='}")
            raise SystemExit()
        else:
            return False


    def path_exists(self,path, v = True):
        '''
        ensure files are present, if so continues if not quits with FileNotFoundError
        input:
            :path: patht to files for pipeline
            :v: if v == True print message, else just check
        output:
            returns True (or fails with FileNotFoundError)
        '''
        
        if not path.exists():
            self.log_messages('warning', f"The {path.name} does not exist.")
            raise FileNotFoundError(f"{path.name}")
        else:
            if v == True:
                self.log_messages('info', f"Found {path.name}.")

            return True

    
        


    def _name_exists(self, name):
        '''
        check if the name is an empty string JOB id can not be empty
       
        '''
        
        if isinstance(name, str):
            if len(name) == 0:
                self.log_messages('warning', 'Job id ca not be empty, please set -j job_id to try again')
                raise SystemExit()
            else:
                return name
        else:
            self.log_messages('warning', 'Job id ca not be empty, please set -j job_id to try again')
            raise SystemExit()

    def link_reads(self, read_source, isolate_id, r_pair):
        '''
        check if read source exists if so check if target exists - if not create isolate dir and link. If already exists report that a dupilcation may have occured in input and procedd

        '''
        # check that job directory exists
        J = pathlib.Path(self.workdir, self.job_id)
        if not J.exists():
            J.mkdir()
        # check that READS exists
        R = J / 'READS'
        if not R.exists():
            R.mkdir()
        
        if f"{read_source}"[0] != '/':
            read_source = self.workdir / read_source
        
        if read_source.exists():
            I = R / f"{isolate_id}" # the directory where reads will be stored for the isolate
            if not I.exists():
                I.mkdir()
            read_target = I / f"{r_pair}"
            if not read_target.exists():
                read_target.symlink_to(read_source)
        else:
            self.log_messages('warning', f"{read_source} does not seem to a valid path. Please check your input and try again.")
            raise SystemExit()

    def unzip_files(self,path, suffix):
        '''
        if a zipped reference is provided try to unzip and then return the unzipped pathname. If unable to unzip then supply message and exit
        input:
            :path: pathname  of file to unzip string
            :unzipped: unzipped path
        '''
        target = self.workdir / path.name.strip(suffix)
        if suffix == '.zip':
            cmd = f"unzip {path} -d {target}"
        elif suffix == 'gz':   
            cmd = f"gzip -d -c {path} > {target}"
        else:
            self.log_messages('warning', f"{path} can not be unzipped. This may be due to file permissions, please provide path to either an uncompressed reference or a file you have permissions to.")
            raise SystemExit

        try:
            subprocess.run(cmd, shell = True)
            return target.name
        except:
            self.log_messages('warning', f"{path} can not be unzipped. This may be due to file permissions, please provide path to either an uncompressed reference or a file you have permissions to.")
            raise SystemExit            

        


    def link_file(self, path):
        '''
        check if file exists and copy to workingdir
        input:
            :path: path to file
        if path does not exist then return false (calling function will request path). 
        if it does exist, then create symlink to the working dir 
        output:
            returns path.name (str)   
        '''
        
            
        if path.exists():
            if f"{path.suffix}" in ['.gz','zip']:
                    self.reference = self.unzip_files(path)
                    if not self.reference.exists():
                        self.log_messages('warning', f"{path} does not exist. Please try again.")
            else:
                target = self.workdir / path.name
                # use rename to copy reference to working directory
                # if the reference is not already in the working directory symlink it to working dir
                if not target.exists():
                    logging.info(f"Linking {path.name} to {self.workdir.name}")
                    target.symlink_to(path)
                    
                    # TODO add in unzip option unzip 
                    found = True
                    
        else:
            self.log_messages('warning', f"Path to {path} does not exist or is not a valid file type (.gbk, .fa, .fasta, .gbk.gz, .fa.gz, .fasta.gz). Please provide a valid path to a file and try again")
            raise SystemExit
            # path = pathlib.Path(path)
        
        return  path.name

    def check_mask(self, mask, original_mask = False):
        '''
        input:
            :mask: path to mask file (str)
            :original_mask: name of orignal mask to use in rerun
        output:
            :mask: path to mask  file in workingdir (str) and a boolean True == rerun snippy-core, tree and distance, False == no need to rerun snippy-core and tree
        '''
        
        # if there is a file path added the generate a symlink
        if len(mask) > 0:
                m = pathlib.Path(mask)
                if f"{m.name}" == original_mask:
                    self.mask = original_mask
                    return original_mask
                else:
                    m = self.link_file(m)
                    self.mask = m
                    return m
        elif len(mask) == 0 and original_mask:
            self.mask = original_mask
            return original_mask
        else:
            return ''
    
    def min_four_samples(self, tab):
        '''
        Ensure that there are a minimum of four samples
        returns True if four or more samples
        '''
        return tab.shape[0] < 4

    def three_cols(self, tab):
        '''
        Ensure that there are 3 columns, isolate, R1 and R2
        returns True if 3 columns False otherwise
        
        '''
        if tab.shape[1] == 3:
            return True
        else:
            return False

    def all_data_filled(self, tab):
        '''
        Ensure that all fields contain data - no NA's
        returns True if there are no nan, False otherwise
        '''
        return tab.isnull().sum().sum() == 0
    

    def check_input_structure(self, tab):
        '''
        check that the structure of the input file is correct, 3 columns, with a minimum of 4 isolates
        :input
            :tab: dataframe of the input file
        :output
            True if make it to the end without exiting with a TypeWarning the user of the reason file was rejected
        '''
        # if the structure of the file is incorrect tell user and kill process
        # not the right information (not three columns)
        
        if not self.three_cols(tab):
            self.log_messages('warning',f"{self.input_file} does not appear to be in the correct configuration")
            raise TypeError(f"{self.input_file} has incorrect number of columns")
        # if there are not enough isolates (>4)
        
        if self.min_four_samples(tab):
            self.log_messages('warning', f"{self.input_file} does not contain enough isolates. The minimum is 4.")
            raise TypeError(f"{self.input_file} has incorrect number of isolates")
        # if any na present indicates that not the full info has been provided
        if not self.all_data_filled(tab):
            self.log_messages('warning',f"{self.input_file} appears to be missing some inforamtion.")
            raise TypeError(f"{self.input_file} appears to be missing some inforamtion.")
        
        return True

    def check_reads_exists(self, tab):
        '''
        check that the correct paths have been given
        if reads not present path_exists will cause a FileNotFound error and warn user
        :input
            :tab: dataframe of the input file
        '''
        for i in tab.itertuples():
            
            if not '#' in i[1]:
                r1 = i[2]
                r2 = i[3]
                self.path_exists(pathlib.Path(r1), v = False)
                self.link_reads(pathlib.Path(r1), isolate_id=f"{i[1].strip()}", r_pair='R1.fq.gz')
                self.path_exists(pathlib.Path(r2), v = False)
                self.link_reads(pathlib.Path(r2), isolate_id=f"{i[1].strip()}", r_pair='R2.fq.gz')
        return True

    def set_isolate_log(self, tab, logfile, validation = False):
        '''
        add the isolates to a log file also adds in anoterh check that this is a rerun of an existing job
        input:
            :tab: dataframe of the isolates to add - same structure as original input file, but not needing to be > 4 isolate
            :logfile: path to logfile
            
        '''        
        self.check_input_structure(tab=tab)
        self.check_reads_exists(tab=tab)
        
        lf = pandas.DataFrame({'Isolate': [i for i in list(tab.iloc[ : , 0]) if '#' not in i ], 'Status': f"INCLUDED", 'Date': self.day})
        lf['Status'] = numpy.where(lf['Isolate'].str.contains('#'), f"REMOVED", lf['Status'])
        isolates = lf[lf['Status'].isin(['INCLUDED', 'ADDED'])]['Isolate']
        lf.to_csv(logfile, sep = '\t', index = False)
        return list(isolates)
        
    
    def set_workflow_input(self, validation = False):
        '''
        read the input file and check that it is the correct format
        input:
            
        output:
            a list of isolates that will be used in generation of job configfile.
        '''
        logfile = pathlib.Path(self.workdir, 'isolates.log') 
        # make df of input file
        

        # the path to the user provided input file check if it exists
        # read the file
        
        tab = pandas.read_csv(self.input_file, sep = None, engine = 'python', header = None)
        
        # self.log_messages('info', f"input file has been opened as a dataframe")
        isolates = self.set_isolate_log(tab = tab, logfile = logfile, validation = validation)
        
        # get a list of isolates
        return(list(set(isolates))) 
        
    def write_pipeline_job(self, maskstring,  script_path = f"{pathlib.Path(__file__).parent / 'utils'}", resource_path = f"{pathlib.Path(__file__).parent / 'templates'}"):
        '''
        write out the pipeline string for transfer to job specific pipeline
        '''
        wd = self.workdir / self.job_id
        
        c = MakeWorkflow()
        config_params = c.write_config_params()
        # config_params = config_params.write_config_params()

        p = MakeWorkflow()
        self.log_messages('message', f"writing the pipeline for {self.job_id}")
        
        # pipeline always has seqdata in 
        pipelinelist = [p.write_all(self.pipeline), p.write_seqdata(), p.write_estimate_coverage(),p.write_generate_yield(script_path),p.write_combine_seqdata()]

        # for just snps contains snippy, snippy-core, dists and tree
        snps_pipeline = [p.write_snippy(), p.write_qc_snippy_initial(), p.write_snippy_core(mask = maskstring), p.write_snp_dists(), p.write_tree(script_path=script_path, alntype='core')]
        # for just assemblies
        assembly_pipeline = [p.write_assemblies(prefillpath = self.prefillpath), p.write_resistome(), p.write_mlst(),p.write_kraken(prefillpath = self.prefillpath), p.write_combine(), p.write_assembly_stats(script_path), p.write_prokka(), p.write_gff_summary(), p.write_combine_kraken()]
        # for roary
        roary_pipeline = [p.write_roary(), p.write_pan_graph(script_path = script_path)]

        # pipeline can always ends with report
        report_pipeline = [p.write_report_collation(pipeline = self.pipeline), p.write_html(pipeline = self.pipeline,workdir = self.workdir, resources = resource_path, job_id = self.job_id, script_path = script_path, assembler = self.assembler_dict[self.assembler])]
        
        if self.pipeline == 's':
            pipelinelist.extend(snps_pipeline)
        elif self.pipeline == 'a':
            pipelinelist.extend(assembly_pipeline)
        elif self.pipeline == 'sa':
            pipelinelist.extend(snps_pipeline)
            pipelinelist.extend(assembly_pipeline)
        elif self.pipeline == 'all':
            pipelinelist.extend(snps_pipeline)
            pipelinelist.extend(assembly_pipeline)
            pipelinelist.extend(roary_pipeline)
        pipelinelist.extend(report_pipeline)


        return(wd, config_params, '\n'.join(pipelinelist))


    def setup_workflow(self, isolates, config_name = 'config.yaml', snake_name = 'Snakefile'):
        '''
        generate job specific snakefile and config.yaml
        input:
            :isolates: a list of isolates that need to be included
        '''

        self.log_messages('info',f"Setting up {self.job_id} specific workflow")
        

        if self.gubbins == True:
            gubbins_string = f"""
        'gubbins.aln', 'gubbins.treefile'
            """
        else:
            gubbins_string = " "
        # make a masking string

        if self.mask != '':
            maskstring = f"--mask {self.workdir / self.mask}"
        else:
            maskstring = ''
                
        # read the config file which is written with jinja2 placeholders (like django template language)
        config_template = jinja2.Template(pathlib.Path(self.resources, 'config_snippy.yaml').read_text())
        config = self.workdir / f"{self.job_id}"/ f"{config_name}"
        
        config.write_text(config_template.render(reference = f"{pathlib.Path(self.workdir, self.ref)}", cpus = self.cpus, name = self.job_id,  minperc = self.minaln,now = self.now, maskstring = maskstring, day = self.day, isolates = ' '.join(isolates)))
        
        self.log_messages('info',f"Config file successfully created")

        snk_template = jinja2.Template(pathlib.Path(self.resources, 'Snakefile_base').read_text())
        snk = self.workdir / snake_name

        # write out custom parts of pipeline
        wd, config_params, pipelinestring = self.write_pipeline_job(maskstring = maskstring)
        snk.write_text(snk_template.render(configfile = f"configfile: '{config_name}'",config_params = config_params, pipeline = pipelinestring, workdir=f"workdir: '{wd}'", gubbins_input = gubbins_string)) 
        
        self.log_messages('info',f"Snakefile successfully created")

 
    def run_workflow(self,snake_name = 'Snakefile'):
        '''
        run snp_detection
        set the current directory to working dir for correct running of pipeline
        if the pipeline wroks, return True else False
        '''
        if self.force:
            force = f"-F"
        else:
            force = f""
        os.chdir(self.workdir)
        if self.dryrun:
            cmd = f"snakemake -np -s {snake_name} 2>&1 | tee -a job.log"
        else:
            cmd = f"snakemake -s {snake_name} -j {self.jobs} {force} 2>&1 | tee -a job.log"
            # cmd = f"snakemake -s {snake_name} --cores {self.cpus} {force} "
        wkf = subprocess.run(cmd, shell = True)
        if wkf.returncode == 0:
            return True
        else:
            return False



    def run_pipeline(self):
        '''
        run pipeline, if workflow runs to completion print out a thank you message.
        '''
        # if -f true force a restart
        if self.force:
            self.force_overwrite()
        # check the pipeline setup 
        self.run_checks()
        # update source data in source.log
        self.set_source_log()
        # open the input file and check it is in the minimal correct format 
        isolates = self.set_workflow_input()
        # setup the workflow files Snakefile and config file
        self.setup_workflow(isolates = isolates)
        # run the workflow
        if self.run_workflow():
            # TODO add in cleanup function to remove snakemkae fluff 
            if not self.dryrun:
                self.log_messages('info', f"Report can be found in {self.job_id}")
                self.log_messages('info', f"Process specific log files can be found in process directories. Job settings can be found in source.log") 
            else:
                if self.force:
                    force = f"-F"
                self.log_messages('info', f"snakemake -j {self.jobs} {force} 2>&1 | tee -a job.log")
            self.log_messages('info', f"Have a nice day. Come back soon.") 
            self.log_messages('info',f"{60 * '='}")

