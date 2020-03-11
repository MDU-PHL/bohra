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
import toml
import shutil
from packaging import version
from Bio import SeqIO, Phylo
from bohra.SnpDetection import RunSnpDetection
# from bohra.bohra_logger import logger
class UpdateBohra(RunSnpDetection):
    '''
    A class to update older bohra directories directories to new bohra structure
    '''
    def __init__(self,args):
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
        self.job_id = self._name_exists(args.job_id)
        self.reference = pathlib.Path(args.reference) if args.reference != '' else ''
        if args.input_file == '':
                self.logger.warning('Input file can not be empty, please set -i path_to_input to try again')
                raise SystemExit()
        else:
            self.input_file = pathlib.Path(args.input_file)
    
    def _update_timestamps(self):
        '''
        if snippy files or assembly present update timestamp so they will not be rerun after linking reads to directory
        '''
        self.logger.info(f"Updating timestamps to prevent uneccesary reruns of snippy and assembler. If you would like to force a rerun use the -F option.")
        for i in self.isolates:
            p = self.jobdir / i
            for f in p.iterdir():
                if "snps" in f"{f}" or f"{f}".endswith(".fa"):
                    cmd = f"touch -m {f}"
                    subprocess.run(cmd, shell = True, capture_output=True)


    def rm_old_reads(self):
        '''
        if READS directory present remove it
        '''

        rds = self.jobdir / 'READS'
        if rds.exists():
            subprocess.run(f"rm -r {rds}", shell = True, capture_output= True)
            # rds.rmdir()
    def archive_old_report(self):
        '''
        if report directory present rename it
        '''

        r = self.jobdir / 'report'
        if r.exists():
            subprocess.run(f"mv {r} {r}_archive", shell = True, capture_output= True)
            # rds.rmdir()
    def _find_reads(self):
        '''
        a function to find sample directories
        '''
        tab = pandas.read_csv(self.input_file, sep = None, engine = 'python', header = None)

        self.jobdir = self.workdir / self.job_id
        self.logger.info(f"Looking for {self.jobdir}")
        if self.jobdir.exists():
            self.logger.info(f"{self.jobdir} found!")
            self.rm_old_reads()
            self.check_input_structure(tab = tab)
            self.check_reads_exists(tab = tab)
            self.isolates = [i.strip('#') for i in tab.iloc[ : , 0]]
            self.archive_old_report()
        else:
            self.logger.warning(f"{self.job_id} does not exist. Are you sure you want to reconfigure this job? Please check help and try again.")
        # pass
    
    def make_tomls(self):
        self.logger.info(f"Checking if tomls need to be made for SNP detection and/or assemblies. This should avoid rerunning these.")
        for i in self.isolates:
            snps = self.jobdir / i / 'snps.aligned.fa'
            if snps.exists():
                d = {i:
                {'snippy':
                    {
                    'reference': f"{self.reference}",
                    'run_snippy': 'Yes',
                    'alignment': f"{i}/snps.aligned.fa",
                    'vcf': f"{i}/snps.vcf",
                    'cmd': f"snippy --outdir {i} --ref {self.reference} --R1 {i}/R1.fq.gz --R2 {i}/R2.fq.gz --force --cpus 8",
                    }
                }
            }
                with open(f"{ self.jobdir / i / 'snippy.toml'}", 'wt') as f:
                    toml.dump(d, f)
            contigs = self.jobdir / i / 'contigs.fa'
            if contigs.exists():
                c = {i:
                {
                    'assembly':
                    {
                        'done': 'Yes'
                    }
                }}
                with open(f"{ self.jobdir / i / 'assembly.toml'}", 'wt') as f:
                    toml.dump(c, f)


    def update(self):

        self._find_reads()
        self._update_timestamps()
        self.make_tomls()


class Nulla2bohra(UpdateBohra):

    '''
    convert a nullarbor directory to a bohra directory
    '''

    def __init__(self,args):
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
        self.reference = pathlib.Path(args.reference) if args.reference != '' else ''
        self.workdir = pathlib.Path(args.workdir)
        self.job_id = self._name_exists(args.job_id)
        if args.input_file == '':
                self.logger.info('Checking your job folder for an input file')
                self.input_file = self.check_for_input_file()
        else:
            self.input_file = pathlib.Path(args.input_file)

    def check_for_input_file(self):

        p = self.workdir / self.job_id

        infile = p / 'input.tab'
        if infile.exists():
            return f"{infile}"
        else:
            self.logger.warning(f"You have not supplied an input file and none can be found in your {self.job_id} directory. Please supply an input file and try again.")
            raise SystemExit

        
class CheckDeps(RunSnpDetection):

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
        self.kraken_db = args.kraken_db 

    def check_sa(self):
        deps = ['shovill', 'skesa', 'spades.py']
        self.check_quals()
        self.check_snippycore()
        self.check_snpdists()
        self.check_kraken2DB()
        # self.check_iqtree()
        # self.check_assembler()
        self.check_assemble_accesories()
        for d in deps:
            self.check_installation(software = d)


    def check(self):
        """
        for bohra check
        """

        self.check_sa()
        self.check_roary()
        # self.logger.info("Well done, you have all dependencies installed correctly!!")
            