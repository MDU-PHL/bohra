import pathlib
import os,getpass
import pandas
import jinja2
import sh
import logging
import filecmp
import datetime
import numpy
import itertools
import subprocess
import re
from packaging import version
from Bio import SeqIO, Phylo
from bohra.SnpDetection import RunSnpDetection
# from bohra.bohra_logger import logger
class ReRunSnpDetection(RunSnpDetection):
    '''
    A class for a Bohra reurn objects - inherits RunSnpDetection
    '''

    def __init__(self, args):

        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        # create file handler which logs even debug messages
        fh = logging.FileHandler('bohra.log')
        fh.setLevel(logging.INFO)
        # create console handler with a higher log level
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        # create formatter and add it to the handlers
        formatter = logging.Formatter('[%(levelname)s:%(asctime)s] %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)
        # add the handlers to the logger
        self.logger.addHandler(ch)
        self.logger.addHandler(fh)
        
        # get date and time
        self.now = datetime.datetime.today().strftime("%d_%m_%y_%H")
        self.day = datetime.datetime.today().strftime("%d_%m_%y")
        # get the working directory
        self.workdir = pathlib.Path(args.workdir)
        # path to pipeline resources
        self.resources = pathlib.Path(args.resources)
        # path to reference => if args.reference is a string, check that it matches existing
        # set force based on args.. this will be set to true if ref is different and/or snippy version 
        self.force = False
        self.assembler = ""
        self.use_singularity = args.use_singularity
        self.singularity_path = args.singularity_path

        self.run_kraken = False
        self.kraken_db = args.kraken_db
        # get original data 
        self.get_source()
        # Reference mask and snippy
        
        if self.pipeline != 'a':
            self.check_reference(new = args.reference)
            # check dependencies
            self.check_for_snippy()
            self.mask = self.check_mask(args.mask, original_mask = self.original_mask)
        elif self.pipeline == 'a':
            self.snippy_version = ''
            self.ref = ''
            self.mask = ''
        # user
        self.user = getpass.getuser()
        # gubbins TODO add back in later!!
        self.gubbins = args.gubbins
        # cluster settings default to command line, 
        self.cluster = args.cluster
        self.json = args.json
        self.queue = args.queue
        # but check to reset if present
        self.get_cluster_reqs()
        # check for cluster settings
        
        
        self.dryrun = args.dry_run
        self.keep = args.keep
        
        self.assembler_dict = {'shovill': 'shovill', 'skesa':'skesa','spades':'spades.py'}
        
        self.cpus = args.cpus
        self.set_snakemake_jobs()

    def get_cluster_reqs(self):
        '''
        check if new cluster configs are being used if not default to stored
        '''
        logger.info(f"Retrieving cluster settings.")
        cluster_log = self.workdir / 'cluster.log'
        if self.cluster == False: #if there is no cluster setting double check if there is an exisitng log 
            if cluster_log.exists(): #reset settings to reflect 
                df = pandas.read_csv(cluster_log, sep = '\t')
                self.json = pathlib.Path(df.loc[df.index[-1], 'cluster_json'])
                self.queue = f"{df.loc[df.index[-1], 'queue']}"
                self.cluster = True
        else:
            self.check_cluster_reqs() # check that settings are appropriate
        
        


    def get_source(self):

        '''
        open the source.log file and extract reference, mask, date, snippy_version
        '''
        logger.info(f"Retrieving settings and software versions.")
        version_pat = re.compile(r'\bv?(?P<major>[0-9]+)\.(?P<minor>[0-9]+)\.(?P<release>[0-9]+)(?:\.(?P<build>[0-9]+))?\b')
        df = pandas.read_csv('source.log', sep = None, engine = 'python')
        
        self.pipeline = df.loc[df.index[-1], 'Pipeline']
        logger.info(f"Previous pipeline was : {self.pipeline}")
        self.use_singularity = df.loc[df.index[-1], 'singularity']
        logger.info(f"Previous --use-singularity was set to : {self.use_singularity}")
        if self.pipeline != 'a':
            self.original_reference = df.loc[df.index[-1], 'Reference']
            logger.info(f"Previous reference was : {self.original_reference}")
            self.original_mask = df.loc[df.index[-1],'Mask'] if isinstance(df.loc[df.index[-1],'Mask'], str) else ''
            logger.info(f"Previous mask was : {self.original_mask}")
            self.original_snippy_version = version_pat.search(df.loc[df.index[-1], 'snippy_version']) if not self.use_singularity else df.loc[df.index[-1], 'snippy_version']
            logger.info(f"Previous snippy_version was : {self.original_snippy_version}")
        if self.pipeline != 's':
            self.assembler = df['Assembler'].unique()[0]
            logger.info(f"Previous assembler used was : {self.assembler}")
        self.orignal_date = df.loc[df.index[-1], 'Date']
        self.input_file = pathlib.Path(f"{df.loc[df.index[-1], 'input_file']}")
        # print(self.input_file)
        self.job_id = df.loc[df.index[-1], 'JobID']
        self.cpus = df.loc[df.index[-1], 'CPUS']
        self.prefillpath = df.loc[df.index[-1], 'prefillpath']
        self.minaln = df.loc[df.index[-1], 'MinAln']
        self.logger.info(f"This is a re-run of an existing job : {self.job_id}. Previous settings will be used.")
        # return reference, mask, snippy_version, date, input_file, pipeline
        

    def check_reference(self, new):
        '''
        Check that reference used in rerun is the same as the reference used in previous run. If not will need to force new SNP detection
        '''
        # check if refs are the same if not set self.ref to new and change to force, else set ref to original
        logger.info(f"Checking reference.")
        if isinstance(new, str) and len(new) > 0 and len(self.original_reference) > 0:
            new_reference = pathlib.Path(f"{new}")
            if f"{new_reference.name}" == self.original_reference:
                self.ref = pathlib.Path(self.original_reference)
            
            else:
                self.ref = self.link_file(path = new_reference)
                self.force = True
                logger.info(f"You have chosen a different reference from the previous run. Snippy will be forced to run again from the beginning.")
        elif isinstance(new, str) and len(new) == 0 and len(self.original_reference) > 0:
            self.ref = pathlib.Path(self.original_reference)
        else:
            logger.warning('There appears to be something wrong with the reference. You will need to run Bohra using the run command.')
            raise SystemExit


    def check_for_snippy(self):
        '''
        Check the version of Snippy, is different will need to force new SNP detection
        '''
        self.check_setup_files()
        if self.use_singularity:
            logger.info(f"You used singulairty containers to run bohra last time, therefore no need to compare snippy versions.")
        else:
            self.current_snippy_version = self.check_deps()
            logger.info(f"Comapring snippy versions.")
            if self.current_snippy_version.group("major", "minor") != self.original_snippy_version.group("major", "minor"):
                self.force = True
                logger.info(f"You are using a different version of Snippy for this re-run, SNP calling will be repeated.")


    def update_source_log(self):
        '''
        update source.log if user changes parameters
        '''
        logger.info(f"Updating {self.job_id} records.")
        df = pandas.read_csv('source.log', sep = None, engine = 'python')
        # if self.pipeline == 'a':
        snippy_v = f'singularity_{self.day}' if self.use_singularity else self.snippy_version
        data =pandas.DataFrame({'JobID':self.job_id, 'Reference':self.ref,'Mask':self.mask, 'Pipeline': self.pipeline, 'CPUS': self.cpus,'MinAln':self.minaln,'Date':self.day, 'User':self.user,'snippy_version':snippy_v ,'input_file':f"{self.input_file}",'prefillpath': self.prefillpath,'Assembler':self.assembler, 'Gubbins':self.gubbins},index=[0])
        df = df.append(data, sort = True)
        df.to_csv('source.log', index=False, sep = '\t')
    
        
    def run_with_gubbins(self):
        '''
        rename core and distance files
        '''
        if self.gubbins:
            self.logger.info(f"You have chosen to run gubbins. Existing core files will be archived and not removed.")
            corefiles = sorted(pathlib.Path(self.workdir, self.job_id).glob('core*'))
            if corefiles:
                for core in corefiles:
                    new = f"{core}".replace('core', 'core_uncorrected')
                    cmd = f"mv {core} {new}"
                    subprocess.run(cmd,shell = True)
            dists = pathlib.Path(self.workdir,self.job_id, 'distances.tab')
            new_dists = f"{dists}".replace('distances', 'distances_uncorrected')
            if dists.exists():
                cmd = f"mv {dists} {new_dists}"
                subprocess.run(cmd,shell = True)
            self.keep = True
    
    def rerun_report(self):
        '''
        Remove report directory from previous run 
        '''
        # os.chdir(self.workdir)
        
        p1 = pathlib.Path(self.workdir, self.job_id, 'report')
        p2 = pathlib.Path(self.workdir, self.job_id, f"report_{self.orignal_date}")
        if self.keep:
            self.logger.info(f"Archiving previous report files")
            cmd = f"mv {p1} {p2}"
            self.remove_core()
        else:
            self.logger.info("Removing previous report files.")
            cmd = f"if [ -d {p1} ];then rm -r {p1}; fi"
        subprocess.run(cmd, shell = True)
    
    def remove_core(self):
        '''
        Need to remove core_isolates.txt to get snakemake to redo snippy core step
        '''
        self.logger.info(f"Removing previous snippy-core output.")
        corefiles = sorted(pathlib.Path(self.workdir, self.job_id).glob('core*'))
        if corefiles:
            for core in corefiles:
                core.unlink()
        
    def check_singularity_directory(self):
        '''
        Check if the singularity directory is empty
        ''' 
        sing_storage_dir = self.workdir / self.job_id / '.snakemake' / 'singularity'
        containers = sorted(sing_storage_dir.glob("*.si*"))
        if len(containers) != 0:
            self.logger.info(f"You already have singularity containers stored. These will be reused.")
            return True
        else:
            self.logger.info(f"There are no singularity containers present. These will be pulled from {self.singularity_path}.")
            return True

    def rerun_checks(self):
        '''
        check if singularity was used last time if so reuse those settings.
        '''

        df = pandas.read_csv('source.log', sep = None, engine = 'python')

        # if 


    def run_pipeline(self):
        '''
        Rerun the pipeline
        '''
        # self.logger.info(f"Previous --use-singularity is still : {self.use_singularity}")
        if self.use_singularity:
            # check that the .singularity directory is present... if not rerun from scratch
            self.check_singularity_directory()
            self.logger.info(f"You have chosen to run bohra with singularity containers. Good luck")

        else:
            # check if the previous run used singularity - if so check for containers.
            self.run_checks()
        # update source
        self.update_source_log()
        self.run_with_gubbins()
        self.rerun_report()
        # self.remove_core()
        
        isolates = self.set_workflow_input()
        # setup the workflow files Snakefile and config file
        self.setup_workflow(isolates = isolates)
        # run the workflow
        if self.run_workflow(): 
            if not self.dryrun:
                self.logger.info(f"Report can be found in {self.job_id}")
                self.logger.info(f"Process specific log files can be found in process directories. Job settings can be found in source.log") 
    
            self.logger.info(f"Have a nice day. Come back soon.") 