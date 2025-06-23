from bohra.launcher.Utils import CustomFormatter, _check_path, _run_subprocess, _get_required_columns, _get_columns_list

import pandas as pd
import pathlib
import os
import logging
import json
import subprocess
import psutil


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




def _get_profile(profile_config: str = '',
                 profile:str = "lcl") -> str:

    # get hostname
    LOGGER.info(f"Tyring to find your profile.")
    # profile = 'lcl'
    if profile_config != '' and _check_path(profile_config):
        with open(profile_config, 'r') as j:
            _cfg = json.load(j)

        p = subprocess.run('hostname', shell = True, capture_output = True, encoding = "utf-8")
        host = p.stdout.strip()
        LOGGER.info(f"Host is : {host}")
        if host in _cfg:
            profile = _cfg[host]
    # if profile == 'no_config':
    #     LOGGER.critical(f"It seems you on an MDU system and trying to run on a head node. Please move to a compute node and try again.")
    #     raise SystemExit
    LOGGER.info(f"You are running bohra with the {profile} profile.")
    return profile
        
def _set_cpu_limit_local(cpus):
# should be used when executor is local
    total_cores = os.cpu_count()
    one,five,fifteen = psutil.getloadavg()
    avail = total_cores - max(one,five,fifteen)

    if int(cpus) == 0:
        return avail
    elif int(cpus) < avail:
        return cpus
    else:
        return avail
