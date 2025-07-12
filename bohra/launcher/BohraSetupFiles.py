from bohra.launcher.Utils import CustomFormatter, _check_path, _run_subprocess, _get_required_columns, _get_columns_list

import pandas as pd
import pathlib
import os
import logging
import shutil

# Logger
LOGGER =logging.getLogger(__name__) 
LOGGER.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(CustomFormatter())
fh = logging.FileHandler('bohra_run.log')
fh.setLevel(logging.DEBUG)
formatter = logging.Formatter('[%(levelname)s:%(asctime)s] %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p') 
fh.setFormatter(formatter)
LOGGER.addHandler(ch) 
LOGGER.addHandler(fh)



def _open_input_file(_file:str) -> pd.DataFrame:
        
    """
    Open input file and return a pandas dataframe
    :param _file: path to input file
    :return: pandas dataframe
    """
    LOGGER.info(f"Opening input file {_file}.")
    if _file == "":
        LOGGER.critical(f"Input file cannot be empty. Please try again.")
        raise SystemExit
    elif not _check_path(_file):
        LOGGER.critical(f"Input file {_file} does not exist or is not accessible.")
        raise SystemExit
    else:
        LOGGER.info(f"Backing up original file {_file} to {pathlib.Path(_file)}.bak.")
        shutil.copy(_file, f"{pathlib.Path(_file)}.bak")
        return pd.read_csv(_file, engine='python', sep = None, dtype = str)

def _check_data_format(df:pd.DataFrame) -> pd.DataFrame:
    """
    Check if the input file is in the correct format
    :param df: pandas dataframe
    :return: pd.DataFrame
    """
    
    columns_req = _get_required_columns()
    columns = _get_columns_list()
    # check if all required columns are present
    for col in columns_req["required"]:
        if col not in df.columns.tolist():
            LOGGER.critical(f"Column {col} is missing from input file.")
            raise SystemExit
    # check if all must have columns are present
    mchk = 0
    for col in columns_req["must_have"]:
        if col in df.columns.tolist():
            mchk = mchk + 1
    if mchk == 0:
        LOGGER.critical(f"Must have at least one of {','.join(columns_req['must_have'])} in your input file.")
        raise SystemExit
    # check if all optional columns are present
    for col in columns:
        if col not in df.columns:
            df[col] = "not_supplied"
    
    for col in df.columns:
        if col not in columns:
            columns.append(col)
    df = df[columns]
    LOGGER.info(f"Saving input file with columns {','.join(columns)} to {pathlib.Path('input_checked.tsv')}.")
    df.to_csv('input_checked.tsv', sep='\t', index=False, header=True)

    return df

def _check_sequence_file(_file:str) -> bool:
    """ Check if the sequence file exists and is accessible 
    :param _file: path to sequence file
    :return: boolean
    """
        
    if _file.exists() and os.access(_file, os.R_OK):
        
        LOGGER.info(f"Found {_file.name}.")
        return True
    else:
        LOGGER.warning(f"The {_file.name} does not exist or is not accessible.")
        return False

def _get_input_types(cols_provided: list, cols_required:list) -> dict:

    """
    Get the input types from the columns provided
    :param cols_provided: list of columns provided in the input file
    :return: list of input types
    """
    input_types = {}
    for c in cols_required:
        if c == "assembly":
            input_types[c] = "contigs.fa"
        elif c == "r1":
            input_types[c] = "R1.fastq.gz"
        elif c == "r2":
            input_types[c] = "R2.fastq.gz"
    return input_types

def _make_workdir(_input:pd.DataFrame, workdir:str) -> bool:

    columns = _get_columns_list()
    if _check_path(workdir):
        
        wd = pathlib.Path(workdir)
        input_table = _open_input_file(_file=_input)
        input_table = _check_data_format(df=input_table)
        
        input_types = _get_input_types(cols_provided=input_table.columns.tolist(), cols_required=columns)
        for row in input_table.iterrows():
            
            for input_type in input_types:
                user_supplied = row[1][input_type]
                LOGGER.info(f"Linking {user_supplied} to {wd / row[1][columns[0]]}/{input_types[input_type]}")
                if _check_sequence_file(pathlib.Path(user_supplied)):
                    try:
                        LOGGER.info(f"Creating directory {row[1][columns[0]]} in {workdir}")
                        target = pathlib.Path(wd / row[1][columns[0]] / input_types[input_type])
                        if not target.exists():
                            pathlib.Path(f"{wd / row[1][columns[0]]}").mkdir(parents=True, exist_ok=True)
                            target.symlink_to(user_supplied)
                        else:
                            LOGGER.warning(f"File {target} already exists. Skipping link creation.")
                    except Exception as e:
                        LOGGER.critical(f"Could not link {user_supplied} to {wd / row[1][columns[0]]}/{input_types[input_type]}. Error: {e}")
                        raise SystemExit
    return "input_checked.tsv"

# class RunSnpDetection(object):
#     '''
#     A class for Bohra pipeline object
#     '''
    
#     def __init__(self, args):
        
#         # user
#         self.user = getpass.getuser()
#         # get date and time
#         self.now = datetime.datetime.today().strftime("%d_%m_%y_%H")
#         self.day = datetime.datetime.today().strftime('%Y-%m-%d')
#         self.script_path = f"{pathlib.Path(__file__).parent}"
#         # get the working directory
        
#         self.assembler = args.assembler
#         self.kraken_db = args.kraken_db
#         self.blast_db = args.blast_db
#         self.data_dir = args.data_dir
#         self.mobsuite_db = args.mobsuite_db if args.mobsuite_db != "" else "no_db"
#         self.workdir = pathlib.Path(args.workdir)        
#         LOGGER.info(f"\033[1mBohra is being run in {self.workdir} by {self.user} on {self.day}.\033[0m")
        
#         # path to pipeline resources
#         self.pipeline = args.pipeline
#         self.preview = True if self.pipeline == 'preview' else False
#         LOGGER.info(f"You are running bohra in {self.pipeline} mode.")
#         self.job_id =args.job_id
#         LOGGER.info(f"Job ID is set {self.job_id}")
#         # path to reference and mask
#         self.ref = args.reference
#         # LOGGER.info(f"The reference is {self.ref}")
#         # self.link_file(self.ref)
#         self.mask = args.mask
#         # path to input file
#         self.reads,self.contigs = self._check_input_files(_input = args.input_file, contigs = args.contigs, pipeline = self.pipeline)

#         self.keep = True if args.keep == 'Y' else False
#         # snippy args
#         self.minmap = args.minmap
#         self.mincov = args.mincov
#         self.basequal = args.basequal
#         self.minqual = args.minqual
#         self.minfrac = args.minfrac
        
#         self.gubbins = args.gubbins
#         self.force = args.force        
#         self.cpus = args.cpus
#         self.config = args.config
#         self.profile_config = os.getenv('BOHRA_PROFILES','')
#         self.profile = self._get_profile() if args.profile == '' else args.profile
#         self.run_kraken = 'false'
#         self.no_phylo = args.no_phylo
#         self.abritamr_args = args.abritamr_args
#         self.spades_args = args.spades_args
#         self.proceed = args.proceed
#         self.use_conda = True if args.no_conda == False else False
#         self.conda_path = args.conda_path
        
        
#     def _check_input_files(self, _input, contigs, pipeline):
        
#         if _input == '':
#             LOGGER.critical(f"You must supply and input file. Please see help and try again.")
#             raise SystemExit

#         # if _input == '':
#         #     if pipeline in ['preview', 'snps','phylogeny','default','full']:
#         #         LOGGER.critical(f"You are trying to run the {pipeline} bohra pipeline - you must supply an input file with paths to reads.")
#         #         raise SystemExit
        
#         return pathlib.Path(_input), contigs


#     def _get_profile(self):

#         # get hostname
#         LOGGER.info(f"Tyring to find your profile.")
#         profile = 'lcl'
#         if self.profile_config != '' and self._path_exists(pathlib.Path(self.profile_config)):
#             with open(self.profile_config, 'r') as j:
#                 _cfg = json.load(j)

#             p = subprocess.run('hostname', shell = True, capture_output = True, encoding = "utf-8")
#             host = p.stdout.strip()
#             LOGGER.info(f"Host is : {host}")
#             if host in _cfg:
#                 profile = _cfg[host]
#         if profile == 'no_config':
#             LOGGER.critical(f"It seems you on an MDU system and trying to run on a head node. Please move to a compute node and try again.")
#             raise SystemExit
#         LOGGER.info(f"You are running bohra with the {profile} profile.")
#         return profile
            
#     def _set_cpu_limit_local(self, cpus):

#         total_cores = os.cpu_count()
#         one,five,fifteen = psutil.getloadavg()
#         avail = total_cores - max(one,five,fifteen)

#         if int(cpus) == 0:
#             return avail
#         elif int(cpus) < avail:
#             return cpus
#         else:
#             return avail

#     def _check_config(self, config):

#         if self._path_exists(pathlib.Path(config)):
#             with open(config, 'r') as c:
#                 data = c.read()
#                 if self.profile in data:
#                     LOGGER.info(f"Profile : {self.profile}.")
#                 else:
#                     LOGGER.warning(f"Profile : {self.profile} is not found in {config}. I hope you know what you are doing!")
#                 return config
#         else:
#             LOGGER.critical(f"You have provided the path to an alternative config file, which does not exist. Please try again.")
#             raise SystemExit

#     def _generate_random_string(self):

#         letters = string.ascii_lowercase
#         result_str = ''.join(random.choice(letters) for i in range(10))

#         return result_str
    
#     def setup_for_rerun(self):
        
#         report_path_orig = self.workdir / 'report'
#         report_path_preview =  self.workdir / f'report_archived_{self.day}_{self._generate_random_string()}'
#         if self.keep and report_path_orig.exists():
#             cmd = f"mv {report_path_orig} {report_path_preview}"
#             LOGGER.info(f"Archiving previous report directory...")
#             self._run_subprocess(cmd  = cmd)
        


#     def _path_exists(self,path, v = True):
#         '''
#         ensure files are present, if so continues if not quits with FileNotFoundError
#         input:
#             :path: patht to files for pipeline
#             :v: if v == True print message, else just check
#         output:
#             returns True (or fails with FileNotFoundError)
#         '''
        
#         if path.exists() and os.access(path, os.R_OK):
#             if v == True:
#                 LOGGER.info(f"Found {path.name}.")
#             return True
#         else:
#             # LOGGER.warning(f"The {path.name} does not exist or is not accessible.")
#             return False

    
#     def _name_exists(self, name):
#         '''
#         check if the name is an empty string JOB id can not be empty
       
#         '''
        
#         if isinstance(name, str):
#             if len(name) == 0:
#                 LOGGER.warning('Job id ca not be empty, please set -j job_id to try again')
#                 raise SystemExit()
#             else:
#                 return name
#         else:
#             LOGGER.warning('Job id ca not be empty, please set -j job_id to try again')
#             raise SystemExit()
   

#     def _check_ref(self, ref, pipeline):
        
#         if pipeline in ['snps','default','full','phylogeny']:
#             if self.ref == '':
#                 LOGGER.critical(f"Reference file must be provided. Please try again.")
#                 raise SystemExit
#             elif not self._path_exists(pathlib.Path(self.ref), v = True):
#                 LOGGER.critical(f"A valid reference file must be provided. Please try again.")
#                 raise SystemExit
            
#             LOGGER.info(f"Reference {self.ref} has been found. Will now copy to running directory.")
#             reference = self._copy_files(_file = self.ref)
#             LOGGER.info(f"Checking if reference is a valid reference file.")
#             p = subprocess.run(f"any2fasta {ref}", shell = True, capture_output = True, encoding = "utf-8")
#             if p.returncode == 0:
#                 LOGGER.info(f"Reference is in a valid format.")
            
#             else:
#                 LOGGER.critical(f"There is something wrong with your reference file. Valid file types are .fasta, .gbk, .fasta.gz, .gbk.gz. Please check your inputs and try again.")
#                 raise SystemExit
#             return reference
#         else:
#             return 'no_ref_req'
        
#     def _check_gzip(self,read):
#         """
#         ensure that gz file is true gz
#         input:
#             read: path to read file
#         """
#         LOGGER.info(f"Checking if file is valid gzip.")
#         if f"{read}".endswith('gz'):
#             p = subprocess.run(f"gzip --test {read}", shell= True, capture_output= True, encoding='utf-8')
#             if p.returncode != 0:
#                 LOGGER.critical(f"Something is wrong with {read}. It is not a valid gzipped file. Please check your inputs and try again.")
#                 raise SystemExit
#             return True
#         else:
#             True
   
    
#     def _check_size_file(self, path):
#         '''
#         check the size of a file
#         '''
#         s = path.stat().st_size
#         return s

#     def _check_kraken2_files(self, k2db):
#         '''
#         ensure that kraken2 DB is not empty
#         '''
#         LOGGER.info(f'Checking that {k2db} is a directory, checking that files are not empty')
#         if pathlib.Path(k2db).is_dir():
#                 LOGGER.info(f'Found {k2db}, checking that files are not empty')
#                 kmerfiles = sorted(pathlib.Path(k2db).glob('*'))
#                 s = []
#                 for k in range(len(kmerfiles)):
#                     s.append(self._check_size_file(pathlib.Path(k2db) / kmerfiles[k]))
#                 if 0 not in s:
#                     self.run_kraken = 'true'
#                     return True
#         else:
#             LOGGER.critical(f'{k2db} is not a directory.')
#             return False
        

#     def _check_kraken2DB(self, checking = False):
#         '''
#         ensure that DB is present and not emtpy
#         '''
        
#         LOGGER.info(f'Searching for kraken2 DB:  $KRAKEN2_DEFAULT_DB')
#         if os.getenv('KRAKEN2_DEFAULT_DB'):
#             try: # check that there is an environment value set for kraken2 db
#                 k2db = pathlib.Path(os.environ["KRAKEN2_DEFAULT_DB"])
#                 LOGGER.info(f"You are using the default kraken2 database at : {os.environ['KRAKEN2_DEFAULT_DB']}")
#                 if self._check_kraken2_files(k2db = k2db):
#                     self.kraken_db = f"{k2db}"
#                     # LOGGER.info(f"You are using the deafult kraken2 database at : {os.environ['KRAKEN2_DEFAULT_DB']}")
#                 LOGGER.info(f"Congratulations your kraken database is present and all files are present.")

#             except KeyError:
#                 LOGGER.critical(f"It seems that your settings for the kraken DB are incorrect, speciation will not be performed.")
#                 self.run_kraken = False # don't run kraken
        
#         else:
            
#             if not checking:
#                 if pathlib.Path(self.kraken_db).exists():
#                     LOGGER.info(f"{self.kraken_db} has been found.")
#                     self.check_kraken2_files(k2db = self.kraken_db)
#                 else:
#                     LOGGER.critical(f"It seems that your settings for the kraken DB are incorrect, speciation will not be performed.")           
#                     self.run_kraken = False
#             else:
#                 LOGGER.warning('You do not have a default kraken2 DB. You have 4 options \n1). bohra run --kraken_db to point to a specific directory \n2). do not perform speciation \n3). set a KRAKEN2_DEFAULT_DB environment variable or 4). use bohra krakendb_download.')
        
#     def _check_mobsuite_db(self):
        
#         if self.mobsuite_db != "no_db":
            
#             if pathlib.Path(self.mobsuite_db).exists() and os.access(self.mobsuite_db, os.W_OK) :
#                 LOGGER.info(f"Your mobsuite database is set up and should be ready to go.")
#             else:
#                 LOGGER.critical(f"Your mobsuite database is not configured properly. Please check your setup and try again.")
#                 raise SystemExit
#         else:
#             LOGGER.info(f"You are using the default mobsuite database in your enviornment. May the force be with you.")

#     def _get_db(self, _type):

#         p = subprocess.run(f"mlst -h", shell = True, capture_output = True, encoding = "utf-8")
#         db = [d for d in p.stdout.strip().split('\n') if _type in d]
#         db = db[0].split('\'')[-2].strip('\'')
#         if db != '':
#             LOGGER.info(f"Found : {db}")
#             return db
#         else:

#             LOGGER.critical(f"There seems to be something wrong with your mlst setup. Please check and try again.")
#             raise SystemExit

#     def _check_mlstdb(self):

#         LOGGER.info(f"Checking mlst setup.")
#         if self.blast_db == "" or self.data_dir == "":
#             LOGGER.warning(f"You do not have mlst databases pre-configured the default DB with your installation of mlst will be used.")
#             # self.blast_db = self._get_db('--blastdb')
#             # self.data_dir = self._get_db('--datadir')
#         elif self._path_exists(pathlib.Path(self.blast_db)) and self._path_exists(pathlib.Path(self.data_dir)):
#             LOGGER.info(f"Your mlst databases have been found.")
#         else:
#             LOGGER.critical(f"There seems to be something wrong with your mlst setup. Please check and try again.")
#             raise SystemExit
        
#         return True


#     def _check_installation(self,software):
#         '''
#         Check that software is installed
#         input:
#             :software: the name of the software - must be command line name 
#         '''

#         if shutil.which(software):
#             LOGGER.info(f"{software} detected.")
#         else:
#             LOGGER.critical(f"{software} is not installed, please check dependencies and try again.")
#             raise SystemExit
    
#     def check_dependencies(self, checking):

                
#         software = {
#             'nextflow':'nextflow -v',
#             'snippy': 'snippy --version 2>&1',
#             'snp-dists': 'snp-dists -v 2>&1',
#             'iqtree': 'iqtree --version 2>&1',
#             'kraken2': 'kraken2 -v 2>&1',
#             'seqkit': 'seqkit version 2>&1',
#             'any2fasta': 'any2fasta -v',
#             'gubbins': 'run_gubbins.py --version 2>&1',
#             'csvtk':'csvtk version 2>&1',
#             'panaroo': 'panaroo --version 2>&1',
#             'shovill': 'shovill -v 2>&1',
#             'skesa': 'skesa -v 2>&1',
#             'spades': 'spades.py -v 2>&1',
#             'abritamr': 'abritamr -v 2>&1',
#             'mlst': 'mlst -v 2>&1',
#             'prokka': 'prokka -v 2>&1',
#             'mash': 'mash --version 2>&1',
#             'quicktree': 'quicktree -v 2>&1',
#             'amrfinderplus': 'amrfinder --version 2>&1',
#             'abritamr': 'abritamr -v 2>&1',
#             'mob_recon': 'mob_recon -V 2>&1'
#         }

#         software_versions = [f"Software"] # a list of versions for output
#         LOGGER.info(f"Checking that dependencies are installed and recording version.")
#         version_pat_3 = re.compile(r'\bv?(?P<major>[0-9]+)\.(?P<minor>[0-9]+)(?:\.(?P<release>[0-9]+)*)?(?:\.(?P<build>[0-9]+)*)?\b')
        
#         for sft in software:
#             p = subprocess.run(software[sft], capture_output=True, encoding = "utf-8", shell = True)
#             p = p.stdout
#             v = version_pat_3.search(p.strip())
#             if v:
#                 v = v.group(0)
#                 software_versions.append(f"{sft} {v}")
#                 LOGGER.info(f"{sft} version {v} detected.")              
#             else:
#                 LOGGER.critical(f"{sft} is not installed.")
#                 raise SystemExit
#         LOGGER.info(f"Now checking kraken2 DB")
#         self._check_kraken2DB(checking = checking)
#         LOGGER.info(f"Checking mobsuite database setup.")
#         self._check_mobsuite_db()
#         software_versions.append(f"kraken2 DB: {self.kraken_db}")
#         if not checking:
#             # print(checking)
#             self._check_mlstdb()
#             software_versions.append(f"mlst blast db: {self.blast_db}")
#             software_versions.append(f"mlst data dir: {self.data_dir}")
#         # write software_versions.tab
        
#             LOGGER.info(f"Writing software_versions.txt")
#             if not pathlib.Path('report').exists():
#                 subprocess.run(f"mkdir -p report", shell = True)
#             pathlib.Path('report','software_versions.txt').write_text('\n'.join(software_versions))
#         LOGGER.info(f"\033[1mCongratulations all dependencies are installed.\033[0m")
#         return True 

#     def _check_shape(self, val, reads = True):
#         if reads and val == 3:
#             return True
#         elif val == 2:
#             return True
#         else:
#             return False
    
#     def _check_gzip(self,read):
#         """
#         ensure that gz file is true gz
#         input:
#             read: path to read file
#         """
#         LOGGER.info(f"Checking if file is valid gzip.")
#         if f"{read}".endswith('gz'):
#             p = subprocess.run(f"gzip --test {read}", shell= True, capture_output= True, encoding='utf-8')
#             if p.returncode != 0:
#                 LOGGER.critical(f"Something is wrong with {read}. It is not a valid gzipped file. Please check your inputs and try again.")
#                 raise SystemExit
#             return True
#         else:
#             True
   
#     def _check_reads(self,reads):

#         if reads != '' and self._path_exists(reads):
            
#             tab = pandas.read_csv(reads, sep = '\t', header = None, dtype = str)
#             if self._check_shape(tab.shape[1]):
#                 LOGGER.info(f"File {reads} is in correct format.")
#             else:
#                 LOGGER.critical(f"{reads} is not in the correct format. Please check your inputs and try again.")
#                 raise SystemExit
#             return tab
#         # elif self.pipeline in ['assemble','amr_typing']:
#         #     return pandas.DataFrame()
#         else:
#             LOGGER.critical(f"Something is wrong with {reads}. Please try again.")
#             raise SystemExit


#     def _check_contigs(self, contigs):

#         if contigs != '' and self._path_exists(pathlib.Path(contigs)):
#             tab = pandas.read_csv(contigs, sep = '\t', header = None, dtype = str)
#             if self._check_shape(tab.shape[1], reads = False):
#                 LOGGER.info(f"File {contigs} is in correct format.")
#                 return True
#             else:
#                 LOGGER.critical(f"{contigs} is not in the correct format. Assemblies will be generated")
#                 return False
#         else:
#             LOGGER.info(f"No valid contigs file has been supplied. Assemblies will be generated.")
#             return False

#     def _check_phylo(self, isolates_list):
#         LOGGER.info(f"Checking if phylo needs to be run.")
#         LOGGER.info(f"Your analysis contains {len(isolates_list)}")
#         if self.pipeline == 'preivew' and len(isolates_list) > 2:
#             LOGGER.info(f"You are running in preview mode, a mash tree will be output")
#         if self.pipeline in ['amr_typing', 'assemble']:
#             LOGGER.info(f"Your are running the {self.pipeline} bohra pipeline. No tree will be generated.")
#             self.no_phylo = True
#         elif len(isolates_list) < 3 and self.pipeline in ['phylogeny', 'default', 'all']:
#             LOGGER.warning(f"You have less than 3 isolates in your dataset, no phylogenetic tree will be generated.")
#             self.no_phylo = True
#         elif not self.no_phylo:
#             LOGGER.info(f"A phylogenetic tree will be generated.")
#         elif self.pipeline == 'all' and len(isolates_list) < 4:
#             LOGGER.warning(f"You can not run panaroo with less than 4 isolates.")
#             self.pipeline = 'default'
        
#         # if self.no_phylo:
#         #     LOGGER.info(f"You have selected no_phylo option so no phylogenetic tree will be generated.")
        
#     def _copy_files(self, _file):

#         p = pathlib.Path(_file)
#         name = p.name
#         if pathlib.Path(self.workdir, name).exists():
#             LOGGER.info(f"The file : {name} already exists in the current directory")
#         else:
#             cmd = f"cp {p} {pathlib.Path(self.workdir, name)}"
#             LOGGER.info(f"Running : {cmd}")
#             subprocess.run(cmd, shell = True)

#         return f"{name}"
    
#     def _link_input_files(self, iso_dir,read, target):

#         if read.exists() and not pathlib.Path(f"{iso_dir}/{target}").exists():
#                     subprocess.run(f"ln -sf {read} {iso_dir}/{target}", shell = True)
#         elif not read.exists():
#             LOGGER.critical(f"{read} is not a valid path. All read files must be valid path. Please try again.")
#             raise SystemExit
    
#     def _check_reads_not_equal(self, r1,r2):

#         if f"{r1}" == f"{r2}":
#             LOGGER.critical(f"It seems that {r1} and {r2} are the same name. You must supply paths to R1 and R2. Please check your inputs and try again.")
#             raise SystemExit

#     def _setup_directory(self, reads):
        
#         isolates_list = []
#         if not reads.empty:
#             LOGGER.info(f"Setting up isolate directories.")
#             for row in reads.iterrows():
#                 if not row[1][0].startswith('#'):
#                     isolates_list.append(row[1][0])
#                     iso_dir = self.workdir/ f"{row[1][0]}" 
#                     if not iso_dir.exists():
#                         subprocess.run(f"mkdir {iso_dir}", shell = True)
#                     # for r in [row[1][1],row[1][2]]:
#                     read1 = pathlib.Path(row[1][1])
#                     read2 =pathlib.Path(row[1][2])
#                     self._check_reads_not_equal(r1 = read1, r2 = read2)
#                     target1 = 'R1.fastq.gz'
#                     target2 = 'R2.fastq.gz'
#                     self._link_input_files(iso_dir = iso_dir, read = read1, target = target1)
#                     self._link_input_files(iso_dir = iso_dir, read = read2, target = target2)       
#             # self._check_phylo(isolates_list = isolates_list)
        
#         elif self.contigs != '' and pathlib.Path(self.contigs).exists():
#             tab = pandas.read_csv(self.contigs, sep = '\t', header = None, names = ['Isolate','Path'])
#             for row in tab.iterrows():
#                 if not row[1][0].startswith('#'):
#                     isolates_list.append(row[1][0])
#                     iso_dir = self.workdir/ f"{row[1][0]}" 
#                     if not iso_dir.exists():
#                         subprocess.run(f"mkdir {iso_dir}", shell = True)
#                     # for r in [row[1][1],row[1][2]]:
#                     contig = pathlib.Path(row[1][1])
#                     target = 'contigs.fa'
#                     self._link_input_files(iso_dir = iso_dir, read = contig, target = target)
        
#         else:
#             LOGGER.critical(f"There seems to be a problem with your input files... not isolates can be extracted. Please check you inputs and try again.")
#             raise SystemExit
#         self._check_phylo(isolates_list = isolates_list)
#         LOGGER.info(f"Updating isolate list.")
#         pathlib.Path(f"isolates.list").write_text('\n'.join(isolates_list))
#         return f"isolates.list"
        
#     def _check_snippy_defaults(self, snippy_arg, val):
        
#         _file = pathlib.Path(self.script_path, 'nextflow.config')

#         if _file.exists():
#             with open(f"{_file}", 'r') as fh:
#                 lines = fh.readlines()
#                 for line in lines:
#                     if snippy_arg in line:
#                         default_val = line.split('=')[-1].strip().strip("\"")
#                         if f"{default_val}" == f"{val}":
#                             # LOGGER.info(f"You are using default setting for {snippy_arg}")
#                             return False
#                         else:
#                             LOGGER.info(f"You have elected to set your own value for {snippy_arg}. Good luck, we hope you know what you are doing ;)")
#                             return True

#     def _generate_snippy_params(self):

#         _dict = {
#                 'mapqual': self.minmap,
#                 'basequal': self.basequal,
#                 'mincov':self.mincov,
#                 'minfrac':self.minfrac,
#                 'minqual':self.minqual
#                 }
#         snippy_opts = []

#         for d in _dict:
#             if self._check_snippy_defaults(snippy_arg=d, val = _dict[d]):
#                 snippy_opts.append(f"--{d} {_dict[d]}")
        
#         return ' '.join(snippy_opts)
    
#     def _generate_spades_args(self):

#         if self.assembler == 'spades':
#             return f"--spades_opt '{self.spades_args}'"
#         else:
#             return ''

#     def _generate_cmd(self, mode, run_kraken, kraken2_db,assembler, mask_string, reference, 
#                         isolates, user, day, contigs, run_iqtree, species, cpus, config, profile, gubbins,blast_db,data_dir, job_id):
        
#         stub = f"nextflow -Dnxf.pool.type=sync run {self.script_path}/main.nf"
#         resume = '' if self.force else "-resume"
#         cpu = f'-executor.cpus={int(cpus)}' if cpus != '' else ''
#         config = f'-c {config}' if config != '' else ''
#         conda = '--enable_conda true' if self.use_conda else ''
#         conda_path = f"--conda_path {self.conda_path}" if self.conda_path != '' else ''
#         snippy_opts = self._generate_snippy_params()
#         spades_opts = self._generate_spades_args()

#         parameters = f"--job_id {job_id} --mode {mode} --run_iqtree {run_iqtree} --run_kraken {run_kraken} --kraken2_db {kraken2_db} --assembler {assembler} --mask_string {mask_string} --reference {reference} --contigs_file {contigs} --species {species if species != '' else 'no_species'} --outdir {self.workdir} --isolates {isolates} --user {user} --day {day} --gubbins {gubbins} --blast_db {blast_db} --data_dir {data_dir} --mobsuite_db {self.mobsuite_db} {conda} {conda_path} {snippy_opts} {spades_opts}"
#         options = f"-with-report bohra_{day}_report.html -with-trace -profile {profile} {resume} {cpu} {config} {'-with-conda' if self.use_conda else ''}"

#         cmd = f"{stub} {parameters} {options}"
#         return cmd
        
#     def _run_subprocess(self, cmd):

#         LOGGER.info(f"Running:\n\033[1m{cmd}\033[0m")
#         p = subprocess.run(cmd, shell = True)
#         return p
    
#     def _run_cmd(self, cmd):
        
#         if self.proceed:
#             return self._run_subprocess(cmd= cmd)
#         else:

#             try:
#                 LOGGER.info(f"Please paste the following command to run the pipeline:\n\033[1m{cmd}\033[0m")
#             except NameError:
#                 LOGGER.critical(f"It seems something has gone wrong with your inputs. \
# Please select a mode to run, choices are 'analysis' or 'finish'")
#                 raise SystemExit

#     def _write_details(self,cmd, reference, mask):

#         col1 = ['Reference file', 'Mask file', 'Working directory', 'User', 'Date', 'Pipeline', 'Nextflow command']
#         col2 = [reference, mask,self.workdir,self.user,self.day,self.pipeline, cmd]
#         df = pandas.DataFrame({'detail':col1,'description':col2})
#         LOGGER.info(f"Saving pipeline details")
#         subprocess.run(f"mkdir -p {self.workdir}/report", shell = True)
#         df.to_csv(f"{self.workdir}/report/details.txt", index =False, sep = '\t')

#     def run_pipeline(self):
#         '''
#         run pipeline, if workflow runs to completion print out a thank you message.
#         '''
#         self.setup_for_rerun()
#         # run checks for inputs
#         if self.use_conda == False:
#             LOGGER.warning(f"You are using a pre-configured conda environment - please note the results may be unexpected.")
#             self.check_dependencies(checking=False)
#         else:
#             LOGGER.info(f"You are running with conda - wise decision!! Will now ensure that kraken DB is configured properly.")
#             self._check_kraken2DB(checking = False)
#             LOGGER.info(f"Checking mobsuite database setup.")
#             self._check_mobsuite_db()
#             LOGGER.info(f"Now looking for MLST setup")
#             self._check_mlstdb()
#         # input files
#         reads = self._check_reads(reads = self.reads)
#         contigs = self._check_contigs(contigs = self.contigs)
        
#         # reference
        
#         reference = self._check_ref(ref = self.ref,pipeline=self.pipeline)
#         # mask
#         if self.mask != '' and not self._path_exists(pathlib.Path(self.mask)):
#             LOGGER.critical(f"{self.mask} is not a valid path please try again.")
#             raise SystemExit
#         elif self.mask != '' and  self._path_exists(pathlib.Path(self.mask, v = False)):
#             self.mask = self._copy_files(_file = self.mask)
#             LOGGER.info(f"Mask file {self.mask} has been provided.")
#         else:
#             LOGGER.info(f"No mask file has been provided.")

#         isolates_list = self._setup_directory(reads = reads) 

#         if contigs:
#             contigs_file = self.contigs
#         else:
#             contigs_file = 'no_contigs'
            
#         run_iqtree = 'false' if self.no_phylo else 'true'
#         mask_string = 'no_mask' if self.mask == '' else self.mask
#         if self.config == '':
#             cpus = self._set_cpu_limit_local(self.cpus)
#             config = ''
#         else:
#             config = self._check_config(self.config)
#             cpus = ''
        
         
#         cmd = self._generate_cmd(mode = self.pipeline, run_kraken = self.run_kraken, kraken2_db = self.kraken_db,
#                         contigs = contigs_file, cpus = cpus, config = config, assembler = self.assembler, 
#                         mask_string = mask_string, reference = reference, run_iqtree = run_iqtree,profile = self.profile,
#                         isolates = isolates_list, day = self.day, user = self.user, 
#                         species = self.abritamr_args, gubbins = self.gubbins, 
#                         blast_db = self.blast_db if self.blast_db != '' else 'no_db', 
#                         data_dir = self.data_dir if self.data_dir != '' else 'no_db', job_id = self.job_id)
#         self._write_details(cmd = cmd, reference = reference, mask = mask_string)
#         self._run_cmd(cmd)

