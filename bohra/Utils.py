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
                if "snps" in f"{f}":
                    cmd = f"touch -m {f}"
                elif f"{f}".endswith(".fa"):
                    cmd = f"touch -m {f}"
                subprocess.run(cmd, shell = True, capture_output=True)


    def rm_old_reads(self):
        '''
        if READS directory present remove it
        '''

        rds = self.jobdir / 'READS'
        if rds.exists():
            rds.rmdir()

    def _find_reads(self):
        '''
        a function to find sample directories
        '''
        tab = pandas.read_csv(self.input_file, sep = None, engine = 'python', header = None)

        self.jobdir = self.workdir / self.job_id
        if self.jobdir.exists():
            self.rm_old_reads()
            self.check_input_structure(tab = tab)
            self.check_reads_exists(tab = tab)
            self.isolates = [i.strip('#') for i in tab.iloc[ : , 0]]
        else:
            self.logger.warning(f"{self.job_id} does not exist. Are you sure you want to reconfigure this job? Please check help and try again.")
        # pass
    
    def update(self):

        self._find_reads()
        self._update_timestamps()


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
        self.workdir = pathlib.Path(args.workdir)
        self.job_id = self._name_exists(args.job_id)
        if args.input_file == '':
                self.logger.warning('Input file can not be empty, please set -i path_to_input to try again')
                raise SystemExit()
        else:
            self.input_file = pathlib.Path(args.input_file)
    