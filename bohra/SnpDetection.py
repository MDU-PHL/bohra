import pathlib
import os, getpass, shutil, re
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
from bohra.write_snakemake import MakeWorkflow
from bohra.write_report import Report


class RunSnpDetection(object):
    '''
    Setup working directory and run snippy
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
        self.version_pat = re.compile(r'\bv?(?P<major>[0-9]+)\.(?P<minor>[0-9]+)\.(?P<release>[0-9]+)(?:\.(?P<build>[0-9]+))?\b')
        self.acc_versions = {}

    def log_messages(self, type, message):
        '''
        Will log messages to the screen with CLEO and also add them to job.log
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
        will force pipeline to run in an existing folder
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
        try:
            snippy = subprocess.run(['snippy', '--version'], stderr=subprocess.PIPE)
            snippy = snippy.stderr.decode().strip()
            self.snippy_version = self.version_pat.search(snippy)
            self.log_messages('info', f"Snippy {snippy} found. Good job!")
            self.acc_versions['snippy'] = f"Snippy {snippy}"
            return(self.version_pat.search(snippy))
        except FileNotFoundError:
            self.log_messages('warning', f"snippy is not installed.")
            raise SystemExit
    
    def check_snippycore(self):
        '''
        check for snippy-core
        '''
        self.check_snippy_versions('snippy-core')

    def check_snpdists(self):
        '''
        check for snp-dists
        '''
        self.check_snippy_versions('snp-dists')



    def check_iqtree(self):
        '''
        check iqtree
        '''
        self.check_snippy_versions('iqtree')

    def check_snippy_versions(self,sft):

        try:
            sft = subprocess.run([software, '--version'], stdout=subprocess.PIPE)
            sft = sft.stdout.decode().strip()
            self.log_messages('info', f"{sft} v.{self.version_pat.search(sft)} found.")
            self.acc_versions[sft] = f"{sft} v.{self.version_pat.search(sft)}"
        except FileNotFoundError:
            self.log_messages('warning', f"{sft} is not installed.")
            raise SystemExit


    def check_assembler(self):
        '''
        check version of assembler
        :software: name of software (str)
        '''
        ret = 0
        assembler_dict = {'shovill': 'shovill', 'skesa':'skesa','spades':'spades.py'}
        try:
            sft = subprocess.run(f"{assembler_dict[self.assembler]} --version", shell = True)
            self.assembler_version = f"{assembler_dict[self.assembler]} v.{self.version_pat.search(sft)}"
            self.log_messages('info', f"{self.assembler_version} found. Good job!")
            self.acc_versions[assembler_dict[self.assembler]] =  f"{self.assembler_version}"
        except subprocess.CalledProcessError as e:
            ret = e.retuncode
        if ret == 0:
            return True
        else:
            return False

    def check_assemble_accesories(self):

        accessories = ['mlst', 'kraken2', 'abricate', 'prokka']
        
        for a in accessories:
            if shutil.which(a):
                self.log_messages('info', f"{a} is installed")
                vers = subprocess.run([a, '--version'], stdout=subprocess.PIPE)
                vers = vers.stdout.decode().strip('\n')
                self.acc_versions[a] =  f"{a} {self.version_pat.search(vers)}"
            else:
                self.log_messages('warning', f"Roary is not installed, please check dependencies and try again.")
                raise SystemExit
    
    def check_roary(self):

        if shutil.which('roary'):
            self.log_messages('info', f"Roary is installed")
        else:
            self.log_messages('warning', f"Roary is not installed, please check dependencies and try again.")
            raise SystemExit


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
        elif self.pipeline != 's':
            if not self.check_assembler():
                self.log_messages('warning', f"The chosen assembler {self.assembler} is not installed. Exiting")
                raise SystemExit
            self.check_assemble_accesories()
        if self.pipeline == 'all':
            self.check_roary()
        

    def check_validation(self, validation_type):

            if validation_type == 'both':
                validate = ['snps', 'clusters']
            else:
                if validation_type in ['snps', 'clusters']:
                    validate = validation_type.split(',')
                else:
                    self.log_messages('warning', f"{validation_type} is not a valid input. Correct options can be taken from {' '.join(['snps', 'clusters', 'both'])}. Please try again")
                    raise SystemExit

            return(validate)
    
    

    def run_checks(self):

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
        set the reference, mask and id for tracking and potential can be updated now - no need for update in rerun.
        input: snippy version if the X.Y previous != X.Y current (if rerun) force to rerun snippy
            
        '''   
        new_df = pandas.DataFrame({'JobID':self.job_id, 'Reference':f"{self.ref}",'Mask':self.mask, 
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
        Determines if snippy.tab exists - as a determination of the job being a rerun of previous job
        input:
            :snippy_tab: file path to snippy_tab
        output:
            False if snippy.tab does not exists or warns user and exits if it does

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
            paths to files for pipeline
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
        
        # check that source is path
        if isinstance(read_source, str):
            read_source = pathlib.Path(read_source)
        # check that source exists
        if read_source.exists():
            I = R / f"{isolate_id}" # the directory where reads will be stored for the isolate
            if not I.exists():
                I.mkdir()
            read_target = I / f"{r_pair}"
            if not read_target.exists():
                read_target.symlink_to(read_source)

    def link_file(self, path, type = 'reference'):
        '''
        check if reference or tree exists and copy to workingdir
        input:
            :workdir: path to working directory for pipeline or current directory
            :ref: path to reference
        if path does not exist then return false (calling function will request path). 
        if it does exist, then create symlink to the working dir 
        output:
            returns path.name (str)   
        '''
        found = False
        while not found: 
            
            if path.exists():
                target = self.workdir / path.name
                # use rename to copy reference to working directory
                # if the reference is not already in the working directory symlink it to working dir
                if not target.exists():
                    logging.info(f"Linking {path.name} to {self.workdir.name}")
                    target.symlink_to(path)
                
                found = True
            else:
                path = self.log_messages('warning', f"Path to {type} does not exist.")
                # path = pathlib.Path(path)
        
        return  path.name

    def check_mask(self, mask, original_mask = False):
        '''
        input:
            :workdir: a path object 
            :mask: path to mask file (str)
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
                    m = self.link_file(m, 'mask')
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
            :filename: name of file the dataframe was built from
            :rerun: if rerun is false, then the limits on 4 isolates is not applied
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
        '''
        for i in tab.itertuples():
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
            :workdir: path to working directory
            :day: day of the run for log (datetime.datetime  as string)
            :filename: name of input file for checks
            :logfile: path to logfile
            :rerun: if True then tab can be smaller than 4 isolates
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
            :workdir: str 
            :job_id: str
            :day: date (str)
            :rerun: boolean (False - for use in check_input_srtucture)
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
        
    def write_pipeline_job(self, maskstring,  script_path = f"{pathlib.Path(__file__).parent}", resource_path = f"{pathlib.Path(__file__).parent / 'templates'}"):
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
        pipelinelist = [p.write_all(self.pipeline), p.write_seqdata(prefillpath = self.prefillpath),p.write_combine_seqdata()]

        # for just snps contains snippy, snippy-core, dists and tree
        snps_pipeline = [p.write_snippy(), p.write_qc_snippy_initial(), p.write_snippy_core(mask = maskstring), p.write_snp_dists(), p.write_tree(script_path=script_path, alntype='core')]
        # for just assemblies
        assembly_pipeline = [p.write_assemblies(prefillpath = self.prefillpath), p.write_resistome(), p.write_mlst(),p.write_kraken(prefillpath = self.prefillpath), p.write_combine(), p.write_assembly_stats(), p.write_prokka(), p.write_gff_summary(), p.write_combine_kraken()]
        # for roary
        roary_pipeline = [p.write_roary(), p.write_pan_graph()]

        # pipeline can always ends with report
        report_pipeline = [p.write_report_collation(pipeline = self.pipeline), p.write_html(pipeline = self.pipeline,workdir = self.workdir, resources = resource_path, job_id = self.job_id, script_path = script_path)]
        
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
            :resources: path to templates (str)
            :workdir: working directory (str) to be used as target directory
            :reference: reference name
            :threads: int
            :id: str - job_id
            :minaln: min align % (int)
            :date: date_time (str)
            :mask: path to mask (str)
            day: date (str)
            :gubbins: Y or N
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
            cmd = f"snakemake -s {snake_name} --cores {self.cpus} {force} 2>&1 | tee -a job.log"
            # cmd = f"snakemake -s {snake_name} --cores {self.cpus} {force} "
        wkf = subprocess.run(cmd, shell = True)
        if wkf.returncode == 0:
            return True
        else:
            return False


    def finish_workflow(self):
        '''
        final message at completion of workflow. If workflow goes to completion print 'thanks for coming message'
        '''
    
        self.log_messages('info', f"Report can be found in {self.job_id}")
        self.log_messages('info', f"Process specific log files can be found in process directories. Job settings can be found in source.log") 
        self.log_messages('info', f"Have a nice day. Come back soon.") 
        self.log_messages('info',f"{60 * '='}")


    def run_pipeline(self):
        
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
            self.log_messages('info', f"Report can be found in {self.job_id}")
            self.log_messages('info', f"Process specific log files can be found in process directories. Job settings can be found in source.log") 
            self.log_messages('info', f"Have a nice day. Come back soon.") 
            self.log_messages('info',f"{60 * '='}")

