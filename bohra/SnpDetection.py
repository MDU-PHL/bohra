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
        self.workdir = pathlib.Path(args.workdir)
        self.check = args.check
        if not self.check:
            LOGGER.info(f"\033[1mBohra is being run in {self.workdir} by {self.user} on {self.day}.\033[0m")
            
            # path to pipeline resources
            self.pipeline = args.pipeline
            self.no_phylo = args.no_phylo
            self.preview = True if self.pipeline == 'preview' else False
            LOGGER.info(f"You are running bohra in {self.pipeline} mode.")
            self.job_id = self._name_exists(args.job_id)
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
            # self.gubbins = args.gubbins 
            # other variables
            # min aln 
            self.minaln = args.minaln
            self.mincov = args.mincov
            self.minqual = args.minqual
            
            self.gubbins = args.gubbins
            self.force = args.force        
            self.cpus = args.cpus
            self.config = args.config
            self.profile = args.profile
            # kraken db settings
            self.run_kraken = False
            self.no_phylo = args.no_phylo
            self.abritamr_args = args.abritamr_args

            
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
    def _set_cpu_limit_local(self, cpus):

        total_cores = os.cpu_count()
        one,five,fifteen = psutil.getloadavg()
        avail = total_cores - max(one,five,fifteen)

        if int(cpus) < avail:
            return cpus
        else:
            return avail

    def _check_config(self, config):

        if self._path_exists(config):
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
            LOGGER.info(f"Archiving preview directory...")
            subprocess.run(cmd, shell = True, encoding = "utf-8", capture_output= True)
        elif self.keep == 'N' and report_path_orig.exists():
            cmd = f"rm {report_path_orig}"
            LOGGER.info(f"Removing previous report files")
            subprocess.run(cmd, shell = True, encoding = "utf-8", capture_output= True)

    # def check_rerun(self):
    #     '''
    #     Check if the job is a rerun of an existing job, if so print message informing user and exit it is considered a rerun if there is a report directory present 

    #     '''
    #     LOGGER.info(f'Checking if job is a rerun of existing job.')
    #     report_path = self.workdir / self.job_id / 'report' / 'index.html'
    #     preview_path = self.workdir / self.job_id / 'report' / 'preview_distances.tab'
    #     # if the path is a string convert to Path
    #     if isinstance(report_path, str):
    #         report_path = pathlib.Path(report_path)
    #     if report_path.exists() and not preview_path.exists():
    #         LOGGER.info(f"This appears to be a rerun of an existing job. Previous result will be removed unless you use --keep.")  
    #         self.rerun_report()
    #         self.remove_core()        
    #         return True
    #     elif preview_path.exists():
    #         self.setup_for_rerun()
    #         return False


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
                    self.run_kraken = True
                    return True
        else:
            LOGGER.critical(f'{k2db} is not a directory.')
            return False
        

    def _check_kraken2DB(self):
        '''
        ensure that DB is present and not emtpy
        '''
        LOGGER.info(f'Searching for kraken2 DB {self.kraken_db}')
        if self.kraken_db == 'KRAKEN2_DEFAULT_DB':
            try: # check that there is an environment value set for kraken2 db
                k2db = pathlib.Path(os.environ["KRAKEN2_DEFAULT_DB"])
                LOGGER.info(f"You are using the deafult kraken2 database at : {os.environ['KRAKEN2_DEFAULT_DB']}")
                if self._check_kraken2_files(k2db = k2db):
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
    
    def check_dependencies(self):

        
        software = {
            'snippy': 'snippy --version 2>&1',
            'snp-dists': 'snp-dists -v 2>&1',
            'iqtree': 'iqtree --version 2>&1',
            'kraken2': 'kraken2 -v 2>&1',
            'seqkit': 'seqkit version 2>&1',
            'any2fasta': 'any2fasta -v',
            'gubbins': 'run_gubbins.py --version 2>&1',
            'csvtk':'csvtk version 2>&1',
            'roary': 'roary -a 2>&1 | tail -n1',
            self.assembler: f"{self.assembler}.py -v 2>&1" if self.assembler == 'spades' else f"{self.assembler} -v 2>&1",
            'abritamr': 'abritamr -v 2>&1',
            'mlst': 'mlst -v 2>&1',
            'prokka': 'prokka -v 2>&1',
            'mash': 'mash --version 2>&1',
            'quicktree': 'quicktree -v 2>&1'
        }

        software_versions = [f"Software"] # a list of versions for output
        LOGGER.info(f"Checking that dependencies are installed and recording version.")
        version_pat_3 = re.compile(r'\bv?(?P<major>[0-9]+)\.(?P<minor>[0-9]+)(\.(?P<release>[0-9]+)(?:\.(?P<build>[0-9]+)*))?\b')
        
        for sft in software:
            try:
                p = subprocess.run(software[sft], capture_output=True, encoding = "utf-8", shell = True)
                p = p.stdout.strip()
                v = version_pat_3.search(p.strip())
                v = v.group(0)
                LOGGER.info(f"{sft} version {v} detected.")
                software_versions.append(f"{sft} {v}")
            except FileNotFoundError:
                LOGGER.critical(f"{sft} is not installed.")
                raise SystemExit
        LOGGER.info(f"Now checking kraken2 DB")
        self._check_kraken2DB()
        software_versions.append(f"kraken2 DB {self.kraken_db}")
        # write software_versions.tab
        if not self.check:
            LOGGER.info(f"Writing software_versions.txt")
            if not pathlib.Path(self.job_id, 'report').exists():
                pathlib.Path(self.job_id, 'report').mkdir()
            pathlib.Path(self.job_id, 'report','software_versions.txt').write_text('\n'.join(software_versions))
        LOGGER.info(f"\033[1mCongratulations all dependencies are installed.\033[0m")
        return True 

    def _check_shape(self, val, reads = True):
        if reads and val == 3:
            return True
        elif val == 2:
            return True
        else:
            return False

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

        if self.pipeline == 'preivew' and len(isolates_list) > 2:
            LOGGER.info(f"You are running in preview mode, a mash tree will be output")
        elif len(isolates_list) < 3 and self.pipeline in ['default', 'all']:
            LOGGER.warning(f"You have less than 3 isolates in your dataset, no phylogenetic tree will be generated.")
            self.no_phylo = True
        elif self.pipeline == 'all' and len(isolates_list) < 4:
            LOGGER.warning(f"You can not run roary with less the 4 isolates.")
            self.pipeline = 'default'
        elif self.no_phylo:
            LOGGER.info(f"You have selected no_phylo option so no phylogenetic tree will be generated.")
        else:
            LOGGER.info(f"A phylogenetic tree will be generated.")

    def _setup_directory(self, reads):
        isolates_list = []
        if not pathlib.Path(self.job_id).exists():
            LOGGER.info(f"Creating job directory")
            subprocess.run(f"mkdir {self.workdir / self.job_id}", shell = True)
        LOGGER.info(f"Setting up isolate directories.")
        for row in reads.iterrows():
            if not row[1][0].startswith('#'):
                isolates_list.append(row[1][0])
                iso_dir = self.workdir/ f"{self.job_id}" /f"{row[1][0]}" 
                if not iso_dir.exists():
                    subprocess.run(f"mkdir {iso_dir}", shell = True)
                for r in [row[1][1],row[1][2]]:
                    read = pathlib.Path(r)
                    target = read.name
                    if read.exists() and not pathlib.Path(f"{iso_dir}/{target}").exists():
                        subprocess.run(f"ln -sf {r} {iso_dir}/{target}", shell = True)
                    elif not read.exists():
                        LOGGER.critical(f"{read} is not a valid path. All read files must be valid path. Please try again.")
                        raise SystemExit
        
        self._check_phylo(isolates_list = isolates_list)
        LOGGER.info(f"Updating isolate list.")
        pathlib.Path(f"isolates.list").write_text('\n'.join(isolates_list))
        return f"isolates.list"
        

    def _generate_cmd(self, min_cov, min_aln,min_qscore, mode, run_kraken, kraken2_db,assembler, mask_string, reference, 
                        outdir, isolates, user, day, contigs, run_iqtree, species, cpus, config, profile, gubbins):
        
        stub = f"nextflow {self.script_path}/main.nf"
        resume = '' if self.force else "-resume"
        cpu = f'-e.cpus={int(cpus)}' if cpus != '' else ''
        config = f'-c {config}' if config != '' else ''
        parameters = f"--min_cov {min_cov} --min_qscore {min_qscore} --min_aln {min_aln} --mode {mode} --run_iqtree {run_iqtree} --run_kraken {run_kraken} --kraken2_db {kraken2_db} --assembler {assembler} \
--mask_string {mask_string} --reference {reference} --contigs_file {contigs} --species {species if species != '' else 'no_species'} --outdir {outdir} --isolates {isolates} --user {user} --day {day} \
--gubbins {gubbins}"
        options = f"-with-report {self.job_id}_seqgen_report.html -with-trace -profile {profile} {resume} {cpu} {config}"

        cmd = f"{stub} {parameters} {options}"
        return cmd
        
        
    def _run_cmd(self, cmd):
        
        try:
            LOGGER.info(f"Please paste the following command to run the pipeline:\n\033[1m{cmd}\033[0m")
        except NameError:
            LOGGER.critical(f"It seems something has gone wrong with your inputs.")
            raise SystemExit


    def run_pipeline(self):
        '''
        run pipeline, if workflow runs to completion print out a thank you message.
        '''
        # run checks for inputs
        self.check_dependencies()
        # input files
        reads = self._check_reads(reads = self.reads)
        contigs = self._check_contigs(contigs = self.contigs)
        
        # reference
        if self.ref == '':
            LOGGER.critical(f"Reference file must be provided. Please try again.")
            raise SystemExit
        elif not self._path_exists(self.ref, v = False):
            LOGGER.critical(f"A valid reference file must be provided. Please try again.")
            raise SystemExit
        else:
            LOGGER.info(f"Reference {self.ref} has been found.")
        # mask
        if self.mask != '' and not self._path_exists(pathlib.Path(self.mask)):
            LOGGER.critical(f"{self.mask} is not a valid path please try again.")
        elif self.mask != '' and  self._path_exists(pathlib.Path(self.mask, v = False)):
            LOGGER.info(f"Mask file {self.mask} has been provided.")
        else:
            LOGGER.info(f"No mask file has been provided.")

        isolates_list = self._setup_directory(reads = reads)

        if contigs:
            contigs_file = self.contigs
        else:
            contigs_file = 'no_contigs'
        run_iqtree = False if self.no_phylo else True

        if self.config == '':
            cpus = self._set_cpu_limit_local(self.cpus)
            config = ''
        else:
            config = self._check_config(self.config)
            cpus = ''
        
        cmd = self._generate_cmd(min_cov = self.mincov, min_aln = self.minaln, min_qscore = self.minqual, mode = self.pipeline, 
                        run_kraken = self.run_kraken, kraken2_db = self.kraken_db,contigs = contigs_file, cpus = cpus, config = config,
                        assembler = self.assembler, mask_string = self.mask, reference = self.ref, run_iqtree = run_iqtree,profile = self.profile,
                        outdir = self.job_id, isolates = isolates_list, day = self.day, user = self.user, species = self.abritamr_args, gubbins = self.gubbins)
        self._run_cmd(cmd)