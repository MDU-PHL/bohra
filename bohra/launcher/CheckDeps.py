
from bohra.launcher.Utils import CustomFormatter

import logging
import subprocess
import pathlib
import os


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

    
def check_dependencies(check:str = "install",
                       envs:str=f"{pathlib.Path(__file__).parent.parent}/environments",
                       force_reinstall:str="false", tool:str="all")->int:
    """
    Install bohra dependencies.
    """
    script_path = f"{pathlib.Path(__file__).parent}"
    LOGGER.info(f"Will now try {check} dependencies. Please be patient this may take some time!!... Maybe get coffee.")
    process = subprocess.Popen(['bash', f"{script_path}/bohra_install.sh", f"{envs}", f"{check}", f"{force_reinstall}", f"{tool}"], stdout=subprocess.PIPE, encoding='utf-8')
    while process.poll() is None:
        l = process.stdout.readline().strip() # This blocks until it receives a newline.
        LOGGER.info(f"{l}")

    if process.returncode != 0:
        # LOGGER.info(f"{process.stderr.read()}")
        LOGGER.error(f"Error {check}ing dependencies.")
        raise SystemError
        # return 1
    else:
        LOGGER.info("Dependencies installed successfully.")
        LOGGER.info("Bohra is ready to go!")
        return 0

def _check_databases(db_install:bool=False)->int:
    """
    Check that required databases are installed.
    """
    script_path = f"{pathlib.Path(__file__).parent}"
    LOGGER.info(f"Will now check databases required for Bohra. Please be patient this may take some time!!... Maybe get coffee.")
    db_check_cmd = 'get' if db_install else 'check'
    process = subprocess.Popen(['bash', f"{script_path}/bohra_databases.sh", f"{db_check_cmd}"],stderr=subprocess.STDOUT, stdout=subprocess.PIPE, encoding='utf-8')
    while process.poll() is None:
        l = process.stdout.readline().strip() # This blocks until it receives a newline.
        LOGGER.info(f"{l}")

    if process.returncode != 0:
        LOGGER.error(f"Error checking databases: {process.stderr}")
        # raise SystemError
        # return 1
    else:
        LOGGER.info("Databases checked successfully.")
        LOGGER.info("Bohra is ready to go!")
        return 0