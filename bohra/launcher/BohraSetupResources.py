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

    LOGGER.info(f"Tyring to find your profile.")
    # profile = 'lcl'
    if _check_path(profile_config):
        try:
            with open(profile_config, 'r') as f:
                _cfg = f.read().strip().split()

            for line in _cfg:
                if profile == line:
                    LOGGER.info(f"Found the profile {profile} in the config file {profile_config}.")
                    return profile
        except Exception as e:
            LOGGER.warning(f"The profile {profile} is not found. Using the default profile 'lcl'.")
    else:
        LOGGER.warning(f"The profile config file {profile_config} does not exist. Using the default profile 'lcl'.")
    profile = 'lcl'
    LOGGER.info(f"You are running bohra with the {profile} profile.")
    return profile
        
def _set_cpu_limit_local(cpus: int = 0) -> int:

    total_cores = os.cpu_count()
    one,five,fifteen = psutil.getloadavg()
    avail = total_cores - max(one,five,fifteen)

    if int(cpus) == 0:
        return avail
    elif int(cpus) < avail:
        return cpus
    else:
        return avail
