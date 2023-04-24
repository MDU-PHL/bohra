import pathlib
import os, getpass, shutil, re, psutil
import pandas
# from bohra.bohra import check_deps
# import sh
import logging
# import filecmp
import datetime
# import numpy
# import itertools
import subprocess
import json
from packaging import version
from bohra.CustomLog import CustomFormatter


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
        
        # user
        self.user = getpass.getuser()
        # get date and time
        self.now = datetime.datetime.today().strftime("%d_%m_%y_%H")
        self.day = datetime.datetime.today().strftime('%Y-%m-%d')
        self.script_path = f"{pathlib.Path(__file__).parent}"
        # get the working directory
        
        self.assembler = args.assembler
        self.kraken_db = args.kraken_db
        self.blast_db = args.blast_db
        self.data_dir = args.data_dir
        self.workdir = pathlib.Path(args.workdir)        
        LOGGER.info(f"\033[1mBohra is being run in {self.workdir} by {self.user} on {self.day}.\033[0m")
        
        # path to pipeline resources
        self.pipeline = args.pipeline
        self.preview = True if self.pipeline == 'preview' else False
        LOGGER.info(f"You are running bohra in {self.pipeline} mode.")
        self.job_id =args.job_id
        LOGGER.info(f"Job ID is set {self.job_id}")
        # path to reference and mask
        self.ref = pathlib.Path(args.reference)
        # LOGGER.info(f"The reference is {self.ref}")
        # self.link_file(self.ref)
        self.mask = args.mask
        # path to input file
        if args.input_file == '':
            LOGGER.critical('`read` file can not be empty, please set -r path_to_input to try again')
            raise SystemExit()
        else:
            self.reads = pathlib.Path(args.input_file)
        self.contigs = args.contigs

        self.keep = True if args.keep == 'Y' else False
        # snippy args
        self.minmap = args.minmap
        self.mincov = args.mincov
        self.basequal = args.basequal
        self.minqual = args.minqual
        self.minfrac = args.minfrac
        
        self.gubbins = args.gubbins
        self.force = args.force        
        self.cpus = args.cpus
        self.config = args.config
        self.profile_config = os.getenv('BOHRA_PROFILES','')
        self.profile = self._get_profile() if args.profile == '' else args.profile
        self.run_kraken = 'false'
        self.no_phylo = args.no_phylo
        self.abritamr_args = args.abritamr_args
        self.proceed = args.proceed
        self.use_conda = True if args.no_conda == False else False
        self.conda_path = args.conda_path
        
        
    
    def _get_profile(self):

        # get hostname
        LOGGER.info(f"Tyring to find your profile.")
        profile = 'lcl'
        if self.profile_config != '' and self._path_exists(pathlib.Path(self.profile_config)):
            with open(self.profile_config, 'r') as j:
                _cfg = json.load(j)

            p = subprocess.run('hostname', shell = True, capture_output = True, encoding = "utf-8")
            host = p.stdout.strip()
            LOGGER.info(f"Host is : {host}")
            if host in _cfg:
                profile = _cfg[host]
        if profile == 'no_config':
            LOGGER.critical(f"It seems you on an MDU system and trying to run on a head node. Please move to a compute node and try again.")
            raise SystemExit
        LOGGER.info(f"You are running bohra with the {profile} profile.")
        return profile
            
    def _set_cpu_limit_local(self, cpus):

        total_cores = os.cpu_count()
        one,five,fifteen = psutil.getloadavg()
        avail = total_cores - max(one,five,fifteen)

        if int(cpus) == 0:
            return avail
        elif int(cpus) < avail:
            return cpus
        else:
            return avail

    def _check_config(self, config):

        if self._path_exists(pathlib.Path(config)):
            with open(config, 'r') as c:
                data = c.read()
                if self.profile in data:
                    LOGGER.info(f"Profile : {self.profile}.")
                else:
                    LOGGER.warning(f"Profile : {self.profile} is not found in {config}. I hope you know what you are doing!")
                return config
        else:
            LOGGER.critical(f"You have provided the path to an alternative config file, which does not exist. Please try again.")
            raise SystemExit

    def _remove_core(self):
        '''
        Need to remove core_isolates.txt to get snakemake to redo snippy core step
        '''
        LOGGER.info(f"Removing previous snippy-core output.")
        corefiles = sorted(pathlib.Path(self.workdir, self.job_id).glob('core*'))
        if corefiles:
            for core in corefiles:
                core.unlink()
    
    def setup_for_rerun(self):
        report_path_orig = self.workdir / self.job_id / 'report'
        report_path_preview =  self.workdir / self.job_id / f'report_{self.day}'
        if self.keep == 'Y' and report_path_orig.exists():
            cmd = f"mv {report_path_orig} {report_path_preview}"
            LOGGER.info(f"Archiving previous report directory...")
            subprocess.run(cmd, shell = True, encoding = "utf-8", capture_output= True)
        elif self.keep == 'N' and report_path_orig.exists():
            cmd = f"rm {report_path_orig}"
            LOGGER.info(f"Removing previous report files")
            subprocess.run(cmd, shell = True, encoding = "utf-8", capture_output= True)


    def _path_exists(self,path, v = True):
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
            return True
        else:
            # LOGGER.warning(f"The {path.name} does not exist or is not accessible.")
            return False

    
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
   

    def _check_ref(self, ref):

        LOGGER.info(f"Checking if reference is a valid reference file.")
        p = subprocess.run(f"any2fasta {ref}", shell = True, capture_output = True, encoding = "utf-8")
        if p.returncode == 0:
            LOGGER.info(f"Reference is in a valid format.")
        else:
            LOGGER.critical(f"There is something wrong with your reference file. Valid file types are .fasta, .gbk, .fasta.gz, .gbk.gz. Please check your inputs and try again.")
            raise SystemExit

    def _check_gzip(self,read):
        """
        ensure that gz file is true gz
        input:
            read: path to read file
        """
        LOGGER.info(f"Checking if file is valid gzip.")
        if f"{read}".endswith('gz'):
            p = subprocess.run(f"gzip --test {read}", shell= True, capture_output= True, encoding='utf-8')
            if p.returncode != 0:
                LOGGER.critical(f"Something is wrong with {read}. It is not a valid gzipped file. Please check your inputs and try again.")
                raise SystemExit
            return True
        else:
            True
   
    
    def _check_size_file(self, path):
        '''
        check the size of a file
        '''
        s = path.stat().st_size
        return s

    def _check_kraken2_files(self, k2db):
        '''
        ensure that kraken2 DB is not empty
        '''
        LOGGER.info(f'Checking that {k2db} is a directory, checking that files are not empty')
        if pathlib.Path(k2db).is_dir():
                LOGGER.info(f'Found {k2db}, checking that files are not empty')
                kmerfiles = sorted(pathlib.Path(k2db).glob('*'))
                s = []
                for k in range(len(kmerfiles)):
                    s.append(self._check_size_file(pathlib.Path(k2db) / kmerfiles[k]))
                if 0 not in s:
                    self.run_kraken = 'true'
                    return True
        else:
            LOGGER.critical(f'{k2db} is not a directory.')
            return False
        

    def _check_kraken2DB(self, checking = False):
        '''
        ensure that DB is present and not emtpy
        '''
        
        LOGGER.info(f'Searching for kraken2 DB:  $KRAKEN2_DEFAULT_DB')
        if os.getenv('KRAKEN2_DEFAULT_DB'):
            try: # check that there is an environment value set for kraken2 db
                k2db = pathlib.Path(os.environ["KRAKEN2_DEFAULT_DB"])
                LOGGER.info(f"You are using the default kraken2 database at : {os.environ['KRAKEN2_DEFAULT_DB']}")
                if self._check_kraken2_files(k2db = k2db):
                    self.kraken_db = f"{k2db}"
                    # LOGGER.info(f"You are using the deafult kraken2 database at : {os.environ['KRAKEN2_DEFAULT_DB']}")
                LOGGER.info(f"Congratulations your kraken database is present and all files are present.")

            except KeyError:
                LOGGER.critical(f"It seems that your settings for the kraken DB are incorrect, speciation will not be performed.")
                self.run_kraken = False # don't run kraken
        
        else:
            
            if not checking:
                if pathlib.Path(self.kraken_db).exists():
                    LOGGER.info(f"{self.kraken_db} has been found.")
                    self.check_kraken2_files(k2db = self.kraken_db)
                else:
                    LOGGER.critical(f"It seems that your settings for the kraken DB are incorrect, speciation will not be performed.")           
                    self.run_kraken = False
            else:
                LOGGER.warning('You do not have a default kraken2 DB. You have 4 options \n1). bohra run --kraken_db to point to a specific directory \n2). do not perform speciation \n3). set a KRAKEN2_DEFAULT_DB environment variable or 4). use bohra krakendb_download.')
        

    def _get_db(self, _type):

        p = subprocess.run(f"mlst -h", shell = True, capture_output = True, encoding = "utf-8")
        db = [d for d in p.stdout.strip().split('\n') if _type in d]
        db = db[0].split('\'')[-2].strip('\'')
        if db != '':
            LOGGER.info(f"Found : {db}")
            return db
        else:

            LOGGER.critical(f"There seems to be something wrong with your mlst setup. Please check and try again.")
            raise SystemExit

    def _check_mlstdb(self):

        LOGGER.info(f"Checking mlst setup.")
        if self.blast_db == "" or self.data_dir == "":
            LOGGER.warning(f"You do not have mlst databases pre-configured the default DB with your installation of mlst will be used.")
            # self.blast_db = self._get_db('--blastdb')
            # self.data_dir = self._get_db('--datadir')
        elif self._path_exists(pathlib.Path(self.blast_db)) and self._path_exists(pathlib.Path(self.data_dir)):
            LOGGER.info(f"Your mlst databases have been found.")
        else:
            LOGGER.critical(f"There seems to be something wrong with your mlst setup. Please check and try again.")
            raise SystemExit
        
        return True


    def _check_installation(self,software):
        '''
        Check that software is installed
        input:
            :software: the name of the software - must be command line name 
        '''

        if shutil.which(software):
            LOGGER.info(f"{software} detected.")
        else:
            LOGGER.critical(f"{software} is not installed, please check dependencies and try again.")
            raise SystemExit
    
    def check_dependencies(self, checking):

                
        software = {
            'nextflow':'nextflow -v',
            'snippy': 'snippy --version 2>&1',
            'snp-dists': 'snp-dists -v 2>&1',
            'iqtree': 'iqtree --version 2>&1',
            'kraken2': 'kraken2 -v 2>&1',
            'seqkit': 'seqkit version 2>&1',
            'any2fasta': 'any2fasta -v',
            'gubbins': 'run_gubbins.py --version 2>&1',
            'csvtk':'csvtk version 2>&1',
            'panaroo': 'panaroo --version 2>&1',
            'shovill': 'shovill -v 2>&1',
            'skesa': 'skesa -v 2>&1',
            'spades': 'spades.py -v 2>&1',
            'abritamr': 'abritamr -v 2>&1',
            'mlst': 'mlst -v 2>&1',
            'prokka': 'prokka -v 2>&1',
            'mash': 'mash --version 2>&1',
            'quicktree': 'quicktree -v 2>&1',
            'amrfinderplus': 'amrfinder --version 2>&1',
            'abritamr': 'abritamr -v 2>&1',
            'mob_recon': 'mob_recon -V 2>&1'
        }

        software_versions = [f"Software"] # a list of versions for output
        LOGGER.info(f"Checking that dependencies are installed and recording version.")
        version_pat_3 = re.compile(r'\bv?(?P<major>[0-9]+)\.(?P<minor>[0-9]+)(?:\.(?P<release>[0-9]+)*)?(?:\.(?P<build>[0-9]+)*)?\b')
        
        for sft in software:
            p = subprocess.run(software[sft], capture_output=True, encoding = "utf-8", shell = True)
            p = p.stdout
            v = version_pat_3.search(p.strip())
            if v:
                v = v.group(0)
                software_versions.append(f"{sft} {v}")
                LOGGER.info(f"{sft} version {v} detected.")              
            else:
                LOGGER.critical(f"{sft} is not installed.")
                raise SystemExit
        LOGGER.info(f"Now checking kraken2 DB")
        self._check_kraken2DB(checking = checking)
        software_versions.append(f"kraken2 DB: {self.kraken_db}")
        if not checking:
            print(checking)
            self._check_mlstdb()
            software_versions.append(f"mlst blast db: {self.blast_db}")
            software_versions.append(f"mlst data dir: {self.data_dir}")
        # write software_versions.tab
        
            LOGGER.info(f"Writing software_versions.txt")
            if not pathlib.Path('report').exists():
                subprocess.run(f"mkdir -p report", shell = True)
            pathlib.Path('report','software_versions.txt').write_text('\n'.join(software_versions))
        LOGGER.info(f"\033[1mCongratulations all dependencies are installed.\033[0m")
        return True 

    def _check_shape(self, val, reads = True):
        if reads and val == 3:
            return True
        elif val == 2:
            return True
        else:
            return False
    
    def _check_gzip(self,read):
        """
        ensure that gz file is true gz
        input:
            read: path to read file
        """
        LOGGER.info(f"Checking if file is valid gzip.")
        if f"{read}".endswith('gz'):
            p = subprocess.run(f"gzip --test {read}", shell= True, capture_output= True, encoding='utf-8')
            if p.returncode != 0:
                LOGGER.critical(f"Something is wrong with {read}. It is not a valid gzipped file. Please check your inputs and try again.")
                raise SystemExit
            return True
        else:
            True
   
    def _check_reads(self,reads):

        if self._path_exists(reads):
            
            tab = pandas.read_csv(reads, sep = '\t', header = None)
            if self._check_shape(tab.shape[1]):
                LOGGER.info(f"File {reads} is in correct format.")
            else:
                LOGGER.critical(f"{reads} is not in the correct format. Please check your inputs and try again.")
                raise SystemExit
            return tab
        else:
            LOGGER.critical(f"Something is wrong with {reads}. Please try again.")
            raise SystemExit


    def _check_contigs(self, contigs):

        if contigs != '' and self._path_exists(pathlib.Path(contigs)):
            tab = pandas.read_csv(contigs, sep = '\t', header = None)
            if self._check_shape(tab.shape[1], reads = False):
                LOGGER.info(f"File {contigs} is in correct format.")
                return True
            else:
                LOGGER.critical(f"{contigs} is not in the correct format. Assemblies will be generated")
                return False
        else:
            LOGGER.info(f"No valid contigs file has been supplied. Assemblies will be generated.")
            return False

    def _check_phylo(self, isolates_list):
        LOGGER.info(f"Checking if phylo needs to be run.")
        LOGGER.info(f"Your analysis contains {len(isolates_list)}")
        if self.pipeline == 'preivew' and len(isolates_list) > 2:
            LOGGER.info(f"You are running in preview mode, a mash tree will be output")
        if len(isolates_list) < 3 and self.pipeline in ['default', 'all']:
            LOGGER.warning(f"You have less than 3 isolates in your dataset, no phylogenetic tree will be generated.")
            self.no_phylo = True
        elif not self.no_phylo:
            LOGGER.info(f"A phylogenetic tree will be generated.")
        if self.pipeline == 'all' and len(isolates_list) < 4:
            LOGGER.warning(f"You can not run roary with less the 4 isolates.")
            self.pipeline = 'default'
        # if self.no_phylo:
        #     LOGGER.info(f"You have selected no_phylo option so no phylogenetic tree will be generated.")
        
    def _copy_files(self, _file):

        p = pathlib.Path(_file)
        name = p.name
        if pathlib.Path(self.workdir, name).exists():
            LOGGER.info(f"The file : {name} already exists in the current directory")
        else:
            cmd = f"cp {p} {pathlib.Path(self.workdir, name)}"
            LOGGER.info(f"Running : {cmd}")
            subprocess.run(cmd, shell = True)

        return f"{name}"
    
    def _link_reads(self, iso_dir,read, target):

        if read.exists() and not pathlib.Path(f"{iso_dir}/{target}").exists():
                    subprocess.run(f"ln -sf {read} {iso_dir}/{target}", shell = True)
        elif not read.exists():
            LOGGER.critical(f"{read} is not a valid path. All read files must be valid path. Please try again.")
            raise SystemExit
        
    def _setup_directory(self, reads):
        
        isolates_list = []
        # if not pathlib.Path(self.job_id).exists():
        #     LOGGER.info(f"Creating job directory")
        #     subprocess.run(f"mkdir {self.workdir}", shell = True)
        LOGGER.info(f"Setting up isolate directories.")
        for row in reads.iterrows():
            if not row[1][0].startswith('#'):
                isolates_list.append(row[1][0])
                iso_dir = self.workdir/ f"{row[1][0]}" 
                if not iso_dir.exists():
                    subprocess.run(f"mkdir {iso_dir}", shell = True)
                # for r in [row[1][1],row[1][2]]:
                read1 = pathlib.Path(row[1][1])
                read2 =pathlib.Path(row[1][2])
                target1 = 'R1.fastq.gz'
                target2 = 'R2.fastq.gz'
                self._link_reads(iso_dir = iso_dir, read = read1, target = target1)
                self._link_reads(iso_dir = iso_dir, read = read2, target = target2)
                
        
        self._check_phylo(isolates_list = isolates_list)
        LOGGER.info(f"Updating isolate list.")
        pathlib.Path(f"isolates.list").write_text('\n'.join(isolates_list))
        return f"isolates.list"
        
    def _check_snippy_defaults(self, snippy_arg, val):
        
        _file = pathlib.Path(self.script_path, 'nextflow.config')

        if _file.exists():
            with open(f"{_file}", 'r') as fh:
                lines = fh.readlines()
                for line in lines:
                    if snippy_arg in line:
                        default_val = line.split('=')[-1].strip().strip("\"")
                        if f"{default_val}" == f"{val}":
                            # LOGGER.info(f"You are using default setting for {snippy_arg}")
                            return False
                        else:
                            LOGGER.info(f"You have elected to set your own value for {snippy_arg}. Good luck, we hope you know what you are doing ;)")
                            return True

    def _generate_snippy_params(self):

        _dict = {
                'mapqual': self.minmap,
                'basequal': self.basequal,
                'mincov':self.mincov,
                'minfrac':self.minfrac,
                'minqual':self.minqual
                }
        snippy_opts = []

        for d in _dict:
            if self._check_snippy_defaults(snippy_arg=d, val = _dict[d]):
                snippy_opts.append(f"--{d} {_dict[d]}")
        
        return ' '.join(snippy_opts)
        
    def _generate_cmd(self, mode, run_kraken, kraken2_db,assembler, mask_string, reference, 
                        isolates, user, day, contigs, run_iqtree, species, cpus, config, profile, gubbins,blast_db,data_dir, job_id):
        
        stub = f"nextflow -Dnxf.pool.type=sync run {self.script_path}/main.nf"
        resume = '' if self.force else "-resume"
        cpu = f'-executor.cpus={int(cpus)}' if cpus != '' else ''
        config = f'-c {config}' if config != '' else ''
        conda = '--enable_conda' if self.use_conda else ''
        snippy_opts = self._generate_snippy_params()
        parameters = f"--job_id {job_id} --mode {mode} --run_iqtree {run_iqtree} --run_kraken {run_kraken} --kraken2_db {kraken2_db} --assembler {assembler} --mask_string {mask_string} --reference {reference} --contigs_file {contigs} --species {species if species != '' else 'no_species'} --outdir {self.workdir} --isolates {isolates} --user {user} --day {day} --gubbins {gubbins} --blast_db {blast_db} --data_dir {data_dir} {conda} {snippy_opts}"
        options = f"-with-report bohra_{day}_report.html -with-trace -profile {profile} {resume} {cpu} {config} {'-with-conda' if self.use_conda else ''}"

        cmd = f"{stub} {parameters} {options}"
        return cmd
        
        
    def _run_cmd(self, cmd):
        
        if self.proceed:
            LOGGER.info(f"Running:\n\033[1m{cmd}\033[0m")
            p = subprocess.run(cmd, shell = True)
            return p
        else:

            try:
                LOGGER.info(f"Please paste the following command to run the pipeline:\n\033[1m{cmd}\033[0m")
            except NameError:
                LOGGER.critical(f"It seems something has gone wrong with your inputs. \
Please select a mode to run, choices are 'analysis' or 'finish'")
                raise SystemExit


    def run_pipeline(self):
        '''
        run pipeline, if workflow runs to completion print out a thank you message.
        '''
        # run checks for inputs
        if self.use_conda == False:
            LOGGER.warning(f"You are using a pre-configured conda environment - please note the results may be unexpected.")
            self.check_dependencies(checking=False)
        else:
            LOGGER.info(f"You are running with conda - wise decision!! Will now ensure that kraken DB is configured properly.")
            self._check_kraken2DB(checking = False)
            LOGGER.info(f"Now looking for MLST setup")
            self._check_mlstdb()
        # input files
        reads = self._check_reads(reads = self.reads)
        contigs = self._check_contigs(contigs = self.contigs)
        
        # reference
        if self.ref == '':
            LOGGER.critical(f"Reference file must be provided. Please try again.")
            raise SystemExit
        elif not self._path_exists(self.ref, v = True):
            LOGGER.critical(f"A valid reference file must be provided. Please try again.")
            raise SystemExit
        else:
            LOGGER.info(f"Reference {self.ref} has been found. Will now copy to running directory.")
            reference = self._copy_files(_file = self.ref)
            self._check_ref(ref = reference)
        # mask
        if self.mask != '' and not self._path_exists(pathlib.Path(self.mask)):
            LOGGER.critical(f"{self.mask} is not a valid path please try again.")
            raise SystemExit
        elif self.mask != '' and  self._path_exists(pathlib.Path(self.mask, v = False)):
            self.mask = self._copy_files(_file = self.mask)
            LOGGER.info(f"Mask file {self.mask} has been provided.")
        else:
            LOGGER.info(f"No mask file has been provided.")

        isolates_list = self._setup_directory(reads = reads)

        if contigs:
            contigs_file = self.contigs
        else:
            contigs_file = 'no_contigs'
            
        run_iqtree = 'false' if self.no_phylo else 'true'
        mask_string = 'no_mask' if self.mask == '' else self.mask
        if self.config == '':
            cpus = self._set_cpu_limit_local(self.cpus)
            config = ''
        else:
            config = self._check_config(self.config)
            cpus = ''
        
         
        cmd = self._generate_cmd(mode = self.pipeline, run_kraken = self.run_kraken, kraken2_db = self.kraken_db,
                        contigs = contigs_file, cpus = cpus, config = config, assembler = self.assembler, 
                        mask_string = mask_string, reference = reference, run_iqtree = run_iqtree,profile = self.profile,
                        isolates = isolates_list, day = self.day, user = self.user, 
                        species = self.abritamr_args, gubbins = self.gubbins, blast_db = self.blast_db if self.blast_db != '' else 'no_db', data_dir = self.data_dir if self.data_dir != '' else 'no_db', job_id = self.job_id)
        self._run_cmd(cmd)

# class CheckBohra(RunSnpDetection):

#     def __init__(self):
        
#         # user
#         self.user = getpass.getuser()
#         # get date and time
#         self.now = datetime.datetime.today().strftime("%d_%m_%y_%H")
#         self.day = datetime.datetime.today().strftime('%Y-%m-%d')
#         self.script_path = f"{pathlib.Path(__file__).parent}"
#         self.check = True
#         self.check_dependencies(checking = self.check)


class SetupInputFiles():

    def __init__(self, args):

        self.read_path = args.read_path
        self.isolate_list = args.isolate_list
        self.now = datetime.datetime.today().strftime("%d_%m_%y_%H")
        self.day = datetime.datetime.today().strftime('%Y-%m-%d')

    def _check_path(self,path):

        """
        Check if path provided exists
        :input - path to where reads/contigs/inputfile are
        :output - boolean 
        """
        LOGGER.info(f"Checking for {path}")
        if pathlib.Path(path).exists():
            LOGGER.info(f"{path} exists.")
            return True
        else:
            LOGGER.critical(f"Path provided : {path} does not exist - please try again.")
            return False
    def _extract_isolates(self, isolate_file):
        
        with open(isolate_file, 'r') as f:
            isolates = f.read().strip().split('\n')
        
        return isolates
    
    def _get_isolate_list(self, _dir,  ext= f"*.f*q.gz"):

        all_data = sorted(_dir.rglob(ext))

        iso_found = set()

        for reads in all_data:
            nme = reads.name.split('_')[0]

            iso_found.add(nme)
        
        LOGGER.info(f"Found {len(list(iso_found))} distinct sample names at path {_dir}")
        return list(iso_found)
    
    

    def _glob_reads(self, _dir, isolates):

        iso_found = self._get_isolate_list(_dir = _dir, ext = "*.f*q.gz")
        
        isolist = []
        if isolates != []:

            for i in list(iso_found):
                for j in isolates:
                    if i in j:
                        isolist.append(i)
        else:
            isolist = iso_found

        lines = []

        for iso in isolist:

            reads = sorted(_dir.rglob(f"*{iso}*.f*q.gz"))
            if len(reads) == 2:
                LOGGER.info(f"Now add reads for {iso}")
                lines.append(f"{iso}\t{reads[0]}\t{reads[1]}")
            elif len(reads) <2:
                LOGGER.warning(f"There do not appear to be 2 reads for {iso}. Skipping")
            else:
                LOGGER.warning(f"There appear to be more than 2 reads available for {iso}. Skipping")
        if lines != []:
            LOGGER.info(f"Saving reads file as isolates.tab")
            pathlib.Path('isolates.tab').write_text('\n'.join(lines))
        else:
            LOGGER.warning(f"It appears that no reads have been found. Please check your input and try again.")
        
        return True

        
    def _glob_data(self, path, isolates =[], data_type = 'reads', ext = "f*q.gz"):

        """
        glob path provided for data
        :input - path provided.
                 data type (reads or contigs - default to reads)
                 ext extension of data type fastq is default
        :output - list of data
        """

        _dir = pathlib.Path(path)

        if data_type == 'reads':

            self._glob_reads(_dir = _dir, isolates= isolates)
        
    
    def find_reads(self):

        if self.read_path != '' and self._check_path(self.read_path):
            if self.isolate_list != '': 
                isolates = self._extract_isolates(isolate_file= self.isolate_list) if self._check_path(self.isolate_list) else []
            else:
                isolates = []
            self._glob_data(path = self.read_path, isolates= isolates)
        else:
            LOGGER.critical(f"There seems to be a problem with your read path. Please provide a valid path to the reads you wish to include")
            raise SystemExit


    


