import pathlib
import os, getpass, shutil, re, psutil
import pandas
import sh
import logging
import filecmp
import datetime
import numpy
import itertools
import subprocess
import json
from packaging import version
from CustomLog import CustomFormatter


# Logger
LOGGER =logging.getLogger(__name__) 
LOGGER.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(CustomFormatter())
fh = logging.FileHandler('bohra.log')
fh.setLevel(logging.DEBUG)
formatter = logging.Formatter('[%(levelname)s:%(asctime)s] %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p') 
fh.setFormatter(formatter)
LOGGER.addHandler(ch) 
LOGGER.addHandler(fh)

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
        self.pipeline = args.pipeline
        self.preview = True if self.pipeline == 'preview' else False
        LOGGER.info(f"You are running bohra in {self.pipeline} mode.")
        self.job_id = self._name_exists(args.job_id)
        LOGGER.info(f"Job ID is set {self.job_id}")
        # path to reference and mask
        self.ref = pathlib.Path(args.reference)
        LOGGER.info(f"The reference is {self.ref}")
        self.link_file(self.ref)
        # 
        # (args.mask)
        if args.mask:
            rerun_core = self.check_mask(args.mask)
        else:
            self.mask = ''
        # path to input file
        if args.reads == '':
            LOGGER.critical('`read` file can not be empty, please set -r path_to_input to try again')
            raise SystemExit()
        else:
            self.reads = pathlib.Path(args.reads)
        self.contigs = args.contigs

        self.keep = True if args.keep == 'Y' else False
        self.check_rerun()
        # self.gubbins = args.gubbins 
        # other variables
        # min aln 
        self.minaln = args.minaln
        self.mincov = args.mincov
        # cluster settings        
        self.user = getpass.getuser()
        
        # self.gubbins = False
        self.force = args.force
        self.dryrun = args.dry_run
        
        self.cpus = args.cpus
        # kraken db settings
        self.kraken_db = args.kraken_db
        self.run_kraken = False
        self.assembler = args.assembler
        self.snippy_version = ''
        self.assembler_dict = {'shovill': 'shovill', 'skesa':'skesa','spades':'spades.py'}
        

    
    
    
    def force_overwrite(self):
        '''
        will force pipeline to run in an existing folder - removes isolate and source logs
        '''
        LOGGER.info(f"You have selected to force overwrite an existing job.")
        isolatelog = self.workdir / f"isolates.log"
        sourcelog = self.workdir / f"source.log"
        # joblog = self.workdir / f"job.log"
        LOGGER.info(f"Removing history.")
        if isolatelog.exists():
            isolatelog.unlink()
        if sourcelog.exists():
            sourcelog.unlink()

        return False

    def check_setup_files(self):
        '''
        check that the working directory, resources directory and the input file exist
        '''
        LOGGER.info(f"Checking that all required input files exist.")
        LOGGER.info(f"Checking that {self.workdir} exists.")
        self.path_exists(self.workdir, v = False) 
        LOGGER.info(f"Checking that {self.resources} exists.")
        self.path_exists(self.resources, v = False)
        LOGGER.info(f"Checking that {self.input_file} exists.")
        self.path_exists(self.input_file)
    
    def check_snippy(self):
        '''
        check for snippy
        '''
        LOGGER.info(f"Checking that snippy is installed and recording version.")
        version_pat = re.compile(r'\bv?(?P<major>[0-9]+)\.(?P<minor>[0-9]+)\.(?P<release>[0-9]+)(?:\.(?P<build>[0-9]+))?\b')
        try:
            snippy = subprocess.run(['snippy --version 2>&1'], capture_output=True, encoding = "utf-8", shell = True)
            snippy = snippy.stdout
            # print(snippy.strip())
            # self.snippy_version = version_pat.search(snippy.strip())
            # LOGGER.info(f"Snippy {self.snippy_version} found. Good job!")
            return(version_pat.search(snippy.strip()))
        except FileNotFoundError:
            LOGGER.info(f"snippy is not installed.")
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
            LOGGER.info(f"{software} is installed")
        else:
            LOGGER.warning(f"{software} is not installed, please check dependencies and try again.")
            raise SystemExit


    def check_assembler(self):
        '''
        check version of assembler
        '''
        ret = 0
        
        self.check_installation(self.assembler_dict[self.assembler])

    def check_abritamr(self):

        if shutil.which('abritamr'):
            LOGGER.info(f"abritamr is installed")
        elif shutil.which('abriTAMR'):
            LOGGER.info(f"abritamr is installed")
        else:
            LOGGER.warning(f"abritamr is not installed, resistome will not be determined. Check https://github.com/MDU-PHL/abritamr for installation instructions.")


    def check_assemble_accesories(self):
        '''
        check the assembly accessories mlst, kraken, abricate and prokka
        '''
        accessories = ['mlst', 'kraken2', 'prokka']
        
        for a in accessories:
            self.check_installation(a)
        self.check_abritamr()
    
    def check_roary(self):
        '''
        check roary is installed
        '''

        self.check_installation('roary')

    def check_size_file(self, path):
        '''
        check the size of a file
        '''
        s = path.stat().st_size
        return s

    def check_kraken2_files(self, k2db):
        '''
        ensure that kraken2 DB is not empty
        '''
        LOGGER.info(f'Checking that {k2db} is a directory, checking that files are not empty')
        print(pathlib.Path(k2db))
        if pathlib.Path(k2db).is_dir():
                LOGGER.info(f'Found {k2db}, checking that files are not empty')
                kmerfiles = sorted(pathlib.Path(k2db).glob('*'))
                s = []
                for k in range(len(kmerfiles)):
                    s.append(self.check_size_file(pathlib.Path(k2db) / kmerfiles[k]))
                if 0 not in s:
                    self.run_kraken = True
                    return True
        else:
            LOGGER.info(f'{k2db} is not a directory.')
            return False
        

    def check_kraken2DB(self):
        '''
        ensure that DB is present and not emtpy
        '''
        LOGGER.info(f'Searching for kraken2 DB {self.kraken_db}')
        if self.kraken_db == 'KRAKEN2_DEFAULT_DB':
            try: # check that there is an environment value set for kraken2 db
                k2db = pathlib.Path(os.environ["KRAKEN2_DEFAULT_DB"])
                LOGGER.info(f"You are using the deafult kraken2 database at : {os.environ['KRAKEN2_DEFAULT_DB']}")
                if self.check_kraken2_files(k2db = k2db):
                    self.kraken_db = f"{k2db}"
                    # LOGGER.info(f"You are using the deafult kraken2 database at : {os.environ['KRAKEN2_DEFAULT_DB']}")

            except KeyError:
                self.run_kraken = False # don't run kraken
        
        else:
            LOGGER.info('You are attempting to use a custom kraken2 DB. This is pretty advanced, good luck!')
            if pathlib.Path(self.kraken_db).exists():
                LOGGER.info(f"{self.kraken_db} has been found.")
                self.check_kraken2_files(k2db = self.kraken_db)
            else:
                LOGGER.warning(f"It seems that your settings for the kraken DB are incorrect.")           
        
        if self.run_kraken:
            LOGGER.info(f"Congratulations your kraken database is present")  
        else:
            LOGGER.warning(f"kraken DB is not installed in the expected path. Speciation will not be performed.")

        
    # def run_with_gubbins(self):
    #     '''
    #     rename core and distance files
    #     '''
    #     if self.gubbins:
    #         LOGGER.info(f"You have chosen to run gubbins. Existing core files will be archived and not removed.")
    #         corefiles = sorted(pathlib.Path(self.workdir, self.job_id).glob('core*'))
    #         if corefiles:
    #             for core in corefiles:
    #                 new = f"{core}".replace('core', 'core_uncorrected')
    #                 cmd = f"mv {core} {new}"
    #                 subprocess.run(cmd,shell = True)
    #         dists = pathlib.Path(self.workdir,self.job_id, 'distances.tab')
    #         new_dists = f"{dists}".replace('distances', 'distances_uncorrected')
    #         if dists.exists():
    #             cmd = f"mv {dists} {new_dists}"
    #             subprocess.run(cmd,shell = True)
    #         self.keep = True
    
    # def rerun_report(self):
    #     '''
    #     Remove report directory from previous run 
    #     '''
    #     # os.chdir(self.workdir)
        
    #     p1 = pathlib.Path(self.workdir, self.job_id, 'report')
    #     p2 = pathlib.Path(self.workdir, self.job_id, f"report_archived_{self.day}")
    #     if self.keep:
    #         LOGGER.info(f"Archiving previous report files")
    #         cmd = f"mv {p1} {p2}"
    #         self.remove_core()
    #     else:
    #         LOGGER.info("Removing previous report files.")
    #         cmd = f"if [ -d {p1} ];then rm -r {p1}; fi"
    #     subprocess.run(cmd, shell = True)
    
    def remove_core(self):
        '''
        Need to remove core_isolates.txt to get snakemake to redo snippy core step
        '''
        LOGGER.info(f"Removing previous snippy-core output.")
        corefiles = sorted(pathlib.Path(self.workdir, self.job_id).glob('core*'))
        if corefiles:
            for core in corefiles:
                core.unlink()
    
    def check_quals(self):

        quals = ['mash', 'seqtk']
        for q in quals:
            self.check_installation(software = q)


    def check_sa(self):
        
        self.check_quals()
        self.check_snippycore()
        self.check_snpdists()
        self.check_kraken2DB()
        self.check_iqtree()
        self.check_assembler()
        self.check_assemble_accesories()

        return(self.check_snippy())

    def check_deps(self):
        '''
        check dependencies Snippy, snippy-core, snp-dists, iqtree
        '''
        # TODO check all software tools used and is there a way to check database last update??
        # TODO check assemblers
        LOGGER.info(f"Checking software dependencies")
        # LOGGER.info(f"{self.pipeline}")
        if self.pipeline in ["sa", "preview"]:
            return(self.check_sa())

        if self.pipeline == "all":
            self.check_roary()
            return(self.check_sa())
            
    
    def run_checks(self):
        '''
        Run checks prior to start - checking all software is installed, if this is a rerun and the input files
        '''
        self.check_setup_files()
        # LOGGER.info(f"{self.ref}")
        # LOGGER.info(f"Checking software dependencies")
        self.snippy_version = self.check_deps()
        # check reference
        if self.ref == '':
            LOGGER.warning(f"You are trying call SNPs, a reference file is required. Please try again using '-r path to reference'")
            raise SystemExit
        else:
            LOGGER.info(f"Starting to link reference")
            self.ref = self.link_file(self.ref)


    def set_source_log(self):
        '''
        set the reference, mask and id for tracking and potential.           
        '''   
        
        # TODO add in options for using singularity containers
        # path if using containers.
        snippy_v = f'singularity_{self.day}' if self.use_singularity else self.snippy_version.group()
        LOGGER.info(f"Snippy : {snippy_v} has been added to the job log file.")
        kraken = self.kraken_db if self.run_kraken else ''
        s = True if self.use_singularity else False
        LOGGER.info(f"Recording your settings for job: {self.job_id}")
        new_df = pandas.DataFrame({'JobID':self.job_id, 'Reference':f"{self.ref}",'Mask':f"{self.mask}", 
                                    'MinAln':self.minaln, 'MinCov': self.mincov, 'Pipeline': self.pipeline, 'CPUS': self.cpus, 'Assembler':self.assembler,
                                    'Date':self.day, 'User':self.user, 'snippy_version':snippy_v, 'input_file':f"{self.input_file}",'prefillpath': self.prefillpath, 'cluster': self.cluster,'singularity': s, 'kraken_db':kraken, 'Gubbins': self.gubbins}, 
                                    index=[0], )
        
        source_path = self.workdir / 'job.log'
        if source_path.exists():
            source_df = pandas.read_csv(source_path, '\t')
            source_df = source_df.append(new_df, sort = True)
        else:
            source_df = new_df

        
        source_df.to_csv(source_path , index=False, sep = '\t')
    
    def setup_for_rerun(self):
        report_path_orig = self.workdir / self.job_id / 'report'
        report_path_preview =  self.workdir / self.job_id / 'report_preview'
        cmd = f"mv {report_path_orig} {report_path_preview}"
        LOGGER.info(f"Archiving preview directory...")
        subprocess.run(cmd, shell = True, encoding = "utf-8", capture_output= True)

    def check_rerun(self):
        '''
        Check if the job is a rerun of an existing job, if so print message informing user and exit it is considered a rerun if there is a report directory present 

        '''
        LOGGER.info(f'Checking if job is a rerun of existing job.')
        report_path = self.workdir / self.job_id / 'report' / 'index.html'
        preview_path = self.workdir / self.job_id / 'report' / 'preview_distances.tab'
        # if the path is a string convert to Path
        if isinstance(report_path, str):
            report_path = pathlib.Path(report_path)
        if report_path.exists() and not preview_path.exists():
            LOGGER.info(f"This appears to be a rerun of an existing job. Previous result will be removed unless you use --keep.")  
            self.rerun_report()
            self.remove_core()        
            return True
        elif preview_path.exists():
            self.setup_for_rerun()
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
        
        if path.exists() and os.access(path, os.R_OK):
            if v == True:
                LOGGER.info(f"Found {path.name}.")
            
        else:
            LOGGER.warning(f"The {path.name} does not exist or is not accessible.")
            raise FileNotFoundError(f"{path.name}")

            return True
    
    def _name_exists(self, name):
        '''
        check if the name is an empty string JOB id can not be empty
       
        '''
        
        if isinstance(name, str):
            if len(name) == 0:
                LOGGER.warning('Job id ca not be empty, please set -j job_id to try again')
                raise SystemExit()
            else:
                return name
        else:
            LOGGER.warning('Job id ca not be empty, please set -j job_id to try again')
            raise SystemExit()

    def link_reads(self, read_source, isolate_id, r_pair):
        '''
        check if read source exists if so check if target exists - if not create isolate dir and link. If already exists report that a dupilcation may have occured in input and procedd

        '''
        # check that job directory exists
        # LOGGER.info(f"Checking that reads are present.")
        R = pathlib.Path(self.workdir, self.job_id)
        if not R.exists():
            R.mkdir()
        # check that READS exists
             
        if f"{read_source}"[0] != '/':
            read_source = self.workdir / read_source
        
        if read_source.exists() and os.access(read_source, os.R_OK):
            I = R / f"{isolate_id}" # the directory where reads will be stored for the isolate
            if not I.exists():
                I.mkdir()
            read_target = I / f"{r_pair}"
            if not read_target.exists():
                read_target.symlink_to(read_source)
        else:
            LOGGER.warning(f"{read_source} does not seem to a valid path. Please check your input and try again.")
            raise SystemExit()

    def unzip_files(self,path, suffix):
        '''
        if a zipped reference is provided try to unzip and then return the unzipped pathname. If unable to unzip then supply message and exit
        input:
            :path: pathname  of file to unzip string
            :unzipped: unzipped path
        '''
        LOGGER.info(f"Checking if reference needs to be unzipped")
        target = self.workdir / pathlib.Path(path).name.strip(suffix)
        
        if suffix == '.zip':
            cmd = f"unzip {path} -d {target}"
        elif suffix == '.gz':   
            cmd = f"gzip -d -c {path} > {target}"
        else:
            LOGGER.warning(f"{path} can not be unzipped. This may be due to file permissions, please provide path to either an uncompressed reference or a file you have permissions to.")
            raise SystemExit

        try:
            LOGGER.info(f"Trying to unzip reference.")
            subprocess.run(cmd, shell = True)
            return target.name
        except:
            LOGGER.warning(f"{path} can not be unzipped. This may be due to file permissions, please provide path to either an uncompressed reference or a file you have permissions to.")
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
        J = self.workdir / self.job_id
        if not J.exists():
            J.mkdir()
        LOGGER.info(f"Getting input files from {path}.") 
        if path.exists():
            if f"{path.suffix}" in ['.gz','zip']:
                    path = pathlib.Path(self.unzip_files(path, f"{path.suffix}"))
                    if not path.exists():
                        LOGGER.warning(f"{path} does not exist. Please try again.")
                        raise SystemExit
            else:
                target = J /  path.name
                # LOGGER.info(f"Target is {target}")
                # use rename to copy reference to working directory
                # if the reference is not already in the working directory symlink it to working dir
                if not target.exists():
                    LOGGER.info(f"Copying {path.name} to {J}")
                    # target.symlink_to(path)
                    subprocess.run(f"cp {path} {target}", shell = True)
                    found = True
                    
        else:
            LOGGER.warning(f"Path to {path} does not exist or is not a valid file type (.gbk, .fa, .fasta, .gbk.gz, .fa.gz, .fasta.gz). Please provide a valid path to a file and try again")
            raise SystemExit
            # path = pathlib.Path(path)
        
        return  path.name

    def check_ref_type(self, ref):

        lst = []
        record = SeqIO.parse(ref, 'genbank')
        for r in record:
            lst.append(r)
        if lst != []:
            return 'genbank'
        else:
            record = SeqIO.parse(ref, 'fasta')
            for r in record:
                lst.append(r)
            if lst != []:
                return 'fasta'
            else:
                LOGGER.info(f"It seems your reference file is not a genbank or fasta file format. Please check your inputs and try again.")
        
    def index_reference(self):

        ref = pathlib.Path(self.job_id , 'ref.fa')
        idx = pathlib.Path(self.job_id , 'ref.fa.fai')
        if f"{pathlib.Path(self.ref).suffix}" in ['.gz','zip']:
                    self.ref = pathlib.Path(self.unzip_files(self.ref, f"{self.ref.suffix}"))
        ref_type = self.check_ref_type(ref = self.ref)
        if ref_type == 'genbank':
            LOGGER.info(f"converting {self.ref} to fasta format.")
            SeqIO.convert(f"{self.ref}", 'genbank', f"{ref}", 'fasta')
        else:
            subprocess.run(f"cp {self.ref} {ref}", shell = True)
        record = SeqIO.parse(f"{ref}", "fasta")
        for r in record:
            r.description = ''
            SeqIO.write(r, f"{ref}", "fasta")
        LOGGER.info(f"Indexing reference.")
        subprocess.run(f"samtools faidx {ref}", shell =True)

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
            LOGGER.info(f"Checking that mask file exists.")
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
        LOGGER.info(f"Checking that there are a minimum of 4 isolates.")
        return tab.shape[0] < 4

    def three_cols(self, tab):
        '''
        Ensure that there are 3 columns, isolate, R1 and R2
        returns True if 3 columns False otherwise
        
        '''
        LOGGER.info(f"Checking that input file is the correct structure.")
        if tab.shape[1] == 3:
            return True
        else:
            return False

    def all_data_filled(self, tab):
        '''
        Ensure that all fields contain data - no NA's
        returns True if there are no nan, False otherwise
        '''
        LOGGER.info("Checking that there is no empty fields in the input file.")
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
            logging.warning(f"{self.input_file} does not appear to be in the correct configuration")
            raise TypeError(f"{self.input_file} has incorrect number of columns")
        # if there are not enough isolates (>4)
        
        if self.min_four_samples(tab):
            LOGGER.warning(f"{self.input_file} does not contain enough isolates. The minimum is 4.")
            raise TypeError(f"{self.input_file} has incorrect number of isolates")
        # if any na present indicates that not the full info has been provided
        if not self.all_data_filled(tab):
            LOGGER.warning('warning',f"{self.input_file} appears to be missing some inforamtion.")
            raise TypeError(f"{self.input_file} appears to be missing some inforamtion.")
        
        return True

    def check_reads_exists(self, tab):
        '''
        check that the correct paths have been given
        if reads not present path_exists will cause a FileNotFound error and warn user
        :input
            :tab: dataframe of the input file
        '''
        LOGGER.info(f"Checking that all the read files exist.")
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
        LOGGER.info(f"Recording the isolates used in job: {self.job_id} on {self.day}")
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
        
        
        isolates = self.set_isolate_log(tab = tab, logfile = logfile, validation = validation)
        
        LOGGER.info(f"This job : {self.job_id} contains {len(list(set(isolates)))}")
        return(list(set(isolates)))         

    def json_setup(self, queue_args):
        '''
        Using the json file provided determine the args to be used in command
        '''
        # print(queue_args)
        LOGGER.info(f"Getting settings from {self.json}")
        try:
            with open(self.json) as f:
                json_file = json.load(f)
            if '__default__' in json_file:
                defs = json_file['__default__']
                arg_list = [i for i in defs]
                arg_cluster = []
                for a in arg_list:
                    if a in queue_args and self.queue == 'sbatch':
                        arg_cluster.append(f"{queue_args[a]} {{cluster.{a}}}")
                    elif a in queue_args and self.queue == 'qsub':
                        string = f"{queue_args[a]} {{cluster.{a}}}" if a not in ['time', 'cpus-per-task', 'mem'] else f"{queue_args[a]}{{cluster.{a}}}"
                        arg_cluster.append(string)
                    else:
                        self.log_messages('warning', f'{a} is not a valid option. Please read docs and try again')
                        raise SystemExit
                return ' '.join(arg_cluster)
        except json.decoder.JSONDecodeError:
            LOGGER.warning(f'There is something wrong with your {self.json} file. Possible reasons for this error are incorrect use of single quotes. Check json format documentation and try again.')


    def cluster_cmd(self):
        snake_name = f"{pathlib.Path(__file__).parent / 'utils'/ 'bohra.smk'}"
        wd = f"{self.workdir / self.job_id}"
        queue_args = ""
        LOGGER.info(f"Setting up cluster settings for {self.job_id} using {self.json}")
        if self.queue == 'sbatch':
            queue_args = {'account':'-A' ,'cpus-per-task':'-c',  'time': '--time', 'partition':'--partition', 'mem':'--mem', 'job':'-J'}
            queue_cmd = f'sbatch'
        elif self.queue == 'qsub':
            queue_args = {'account':'-P' ,'cpus-per-task': '-l ncpus=',  'time': '-l walltime=', 'partition':'-q', 'mem':'-l mem=', 'job':'-N'}
            queue_cmd = f'qsub'
        else:
            logging.warning(f'{self.queue} is not supported please select sbatch or qsub. Alternatively contact developer for further advice.')
            raise SystemExit
    
        queue_string = self.json_setup(queue_args = queue_args)

        return f"snakemake -j 999 --cluster-config {self.json} --cluster '{queue_cmd} {queue_string}'"

    
    
    def setup_workflow(self, isolates, config_name = 'config.yaml'):
        '''
        generate job specific snakefile and config.yaml
        input:
            :isolates: a list of isolates that need to be included
        '''

        LOGGER.info(f"Setting up {self.job_id} specific workflow")
        
        # make a masking string
        wd = self.workdir / self.job_id
        if self.mask != '':
            maskstring = f"{self.workdir / self.job_id /self.mask}"
        else:
            maskstring = ''
        LOGGER.info(f'Mask string : {maskstring}')
        LOGGER.info(f"Writing config file for job : {self.job_id}")
        reference = pathlib.Path(self.ref)
        vars_for_file = {
            'workdir': f"{wd}",
            'prefill_path' : self.prefillpath,
            'job_id' : self.job_id,
            'assembler' : self.assembler if self.pipeline != 's' else 'no_assembler',
            'mask_string': maskstring, 
            'template_path':f"{pathlib.Path(__file__).parent / 'templates'}",
            'script_path':f"{pathlib.Path(__file__).parent / 'utils'}",
            'reference' : f"{wd / reference.name}",
            'minperc' : self.minaln,
            'now' : self.now,
            'day': self.day, 
            'isolates' : ' '.join(isolates),
            'gubbins': self.gubbins,
            'pipeline': self.pipeline,
            'min_cov': self.mincov,
            'run_kraken': 'speciation' if self.run_kraken else 'no_speciation',
            'kraken_db': f"{self.kraken_db}",
            'kraken_threads': self.kraken_threads,
            'preview': self.preview, 
            'prefill_path': self.prefillpath if self.prefillpath != '' else 'nopath',
            'snippy_singularity': self.snippy_singularity,
            'abritamr_singularity': self.abritamr_singularity
        }

        # read the config file which is written with jinja2 placeholders (like django template language)
        config_template = jinja2.Template(pathlib.Path(self.resources, 'config_snippy.yaml').read_text())
        config = self.workdir / f"{self.job_id}"/ f"{config_name}"
        
        config.write_text(config_template.render(vars_for_file))
        
        LOGGER.info(f"Config file successfully created")
       
    def run_workflow(self):
        '''
        run snp_detection
        set the current directory to working dir for correct running of pipeline
        if the pipeline works, return True else False
        '''
        snake_name = f"{pathlib.Path(__file__).parent / 'utils'/ 'bohra.smk'}"
        if self.use_singularity:
            singularity_string = f"--use-singularity --singularity-args \"--home /home\""
        else:
            singularity_string = ''

        if self.force:
            force = f"-F"
        else:
            force = f""
        os.chdir(self.workdir)
        
        if self.dryrun:
            dry = '-np'
        else:
            dry = ''
        wd = self.workdir / self.job_id
        if self.cluster:
            cmd = f"{self.cluster_cmd()} -s {snake_name} -d {wd} {force} {singularity_string} --latency-wait 1200"
        else:
            cmd = f"snakemake {dry} -s {snake_name} {singularity_string} -j {self.cpus} -d {self.job_id} {force} {self.verbose} 2>&1 | tee -a bohra.log"
            # cmd = f"snakemake -s {snake_name} --cores {self.cpus} {force} "
        LOGGER.info(f"Running job : {self.job_id} with {cmd} this may take some time. We appreciate your patience.")
        wkf = subprocess.run(cmd, shell = True)
        
        while True:
            if wkf.stdout != None:
                line = wkf.stdout.readline().strip()
                if not line:
                    break
            line = ''
            break
            self.LOGGER.info(f"{line}")
        if wkf.returncode == 0:
            return True
        else:
            return False



    def run_pipeline(self):
        '''
        run pipeline, if workflow runs to completion print out a thank you message.
        '''
        # need to keep track of reference file - use hash?
        # job_id,reference_file,reference_hash
        
        
        # if -f true force a restart
        if self.force:
            self.force_overwrite()
        # check the pipeline setup 
        