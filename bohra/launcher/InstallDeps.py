
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
fh = logging.FileHandler('bohra_installation.log')
fh.setLevel(logging.DEBUG)
formatter = logging.Formatter('[%(levelname)s:%(asctime)s] %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p') 
fh.setFormatter(formatter)
LOGGER.addHandler(ch) 
LOGGER.addHandler(fh)

    
def install_dependencies(prefix):
    """
    Install bohra dependencies.
    """
    script_path = f"{pathlib.Path(__file__).parent}"
    LOGGER.info(f"Will now try to install dependencies. Please be patient this may take some time!!... Maybe get coffee.")
    process = subprocess.Popen(['bash', f"{script_path}/bohra_install.sh", f"{prefix}"], stdout=subprocess.PIPE, encoding='utf-8')
    while process.poll() is None:
        l = process.stdout.readline().strip() # This blocks until it receives a newline.
        print(f"{l}")

    if process.returncode != 0:
        LOGGER.error(f"Error installing dependencies: {process.returncode}")
        raise SystemError
    else:
        LOGGER.info("Dependencies installed successfully.")
        LOGGER.info("Bohra is ready to go!")
        return True
        
        