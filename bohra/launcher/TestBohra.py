import pathlib
import datetime
import logging
from bohra.launcher.Utils import CustomFormatter, _check_path, _run_subprocess
from bohra.launcher.SetupInput import find_data
# Logger
LOGGER =logging.getLogger(__name__) 
LOGGER.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(CustomFormatter())
fh = logging.FileHandler('bohra_test.log')
fh.setLevel(logging.DEBUG)
formatter = logging.Formatter('[%(levelname)s:%(asctime)s] %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p') 
fh.setFormatter(formatter)
LOGGER.addHandler(ch) 
LOGGER.addHandler(fh)


def _download_reads_from_github(isolate_list:list, download_stub:str):
    """
    Downloads test data from github.
    """
    LOGGER.info(f"Downloading test data from github.")

    for isolate in isolate_list:
        LOGGER.info(f"Downloading reads for {isolate}.")
        for r in [1,2]:
            cmd = f"mkdir -p test_data/{isolate} && wget -O test_data/{isolate}/{isolate}_{r}.fastq.gz {download_stub}/{isolate}/{isolate}_{r}.fastq.gz"
            LOGGER.info(f"Running : {cmd}")
            _run_subprocess(cmd = cmd)
        
def _download_reference_from_github(download_stub:str):

    cmd = f"wget -O Lm_Cluster1_J1-108.fa {download_stub}/Lm_Cluster1_J1-108.fa"
    LOGGER.info(f"Downloading Lm_Cluster1_J1-108.fa to {pathlib.Path.cwd() / 'Lm_Cluster1_J1-108.fa'}")
    LOGGER.info(f"Running : {cmd}")
    _run_subprocess(cmd = cmd)
    return 'Lm_Cluster1_J1-108.fa'

def _check_reference_test(path, download_stub):

    if not pathlib.Path(path).exists():

        return _download_reference_from_github(download_stub=download_stub)

    else:
        return path

def _check_test_data(path, isolate_list):

    for i in isolate_list:
        reads = sorted(pathlib.Path(path).glob(f"{i}/{i}*.f*q.gz"))
        if len(reads) != 2:
            LOGGER.warning(f"Reads for {i} are not found. Will try to get them for you!")
            return False
    LOGGER.info(f"All reads are found at {path}")
    return True

def run_tests(cpus:int=1):
    download_stub = "https://raw.githubusercontent.com/MDU-PHL/bohra/master/data"
    read_path = f"{pathlib.Path.cwd() / 'test_data'}"
    isolate_list = ['ERR1102348','ERR1102353','ERR1102355','ERR1102356']
    reference = _check_reference_test(path=f"{pathlib.Path.cwd() / 'Lm_Cluster1_J1-108.fa'}", download_stub=download_stub)
    LOGGER.info(f"Checking availability of data.")
    if not pathlib.Path(read_path).exists():
        LOGGER.info(f"Will now download some reads for testing - this may take a little while - it might be coffee time.")
        _download_reads_from_github(download_stub=download_stub, isolate_list=isolate_list)
        LOGGER.info(f"Reads have been downloaded to {read_path}.")
        LOGGER.info(f"Now generating the input file from the reads.")
        find_data(reads = f"{read_path}",contigs="",isolate_ids ="", outname="bohra_input.tsv" )
        
    elif _check_test_data(path = read_path, isolate_list = isolate_list):
        find_data(reads = f"{read_path}",contigs="",isolate_ids ="", outname="bohra_input.tsv" )
        
    cmd = f"bohra run full -i bohra_input.tsv -ref {reference} --cpus {cpus} --report_outdir bohra_test_output_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
    LOGGER.info(f"Now testing that the bohra installation has worked. Running command: {cmd}")
    proc = _run_subprocess(cmd=cmd)

    if proc.returncode == 0:
        LOGGER.info(f"bohra test has completed successfully!!")
    else:
        LOGGER.critical(f"bohra run was not successful... The following error was reported : {proc.stderr}. Please raise an issue on github.")

