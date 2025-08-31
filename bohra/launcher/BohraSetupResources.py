from bohra.launcher.Utils import CustomFormatter, _check_path, _run_subprocess, _get_required_columns, _get_columns_list

import pandas as pd
import pathlib
import os
import logging
import json
import subprocess
import psutil

from bohra.launcher.Utils import _get_template, _get_target


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



def _generate_nxf_config(title: str, res: dict, wd: str) -> pathlib.Path:
    try:
        template  = f"{pathlib.Path(__file__).parent.parent.resolve() / 'templates' / 'nf.config.j2'}"
        if not _check_path(template):
            LOGGER.critical(f"Template file {template} does not exist. Cannot generate Nextflow config file.")
            raise SystemExit
        

        LOGGER.info(f"Loading template {template} for config generation.")
        template = _get_template(template)
        target = _get_target(wd, title)
    #     print(data["title"])
        print("Rendering config file...")
        target.write_text(template.render(res))
    except Exception as e:
        LOGGER.critical(f"Error generating Nextflow config file: {e}")
        raise SystemExit
    LOGGER.info(f"Nextflow config file generated at {target}.")
    return target


def _get_cpu_limit_local(cpus: int = 0) -> int:

    total_cores = os.cpu_count()
    one,five,fifteen = psutil.getloadavg()
    avail = total_cores - max(one,five,fifteen)

    return avail
    
def _get_default_profile_cfg() -> dict:

    cfg_file = f"{pathlib.Path(__file__).parent.parent.resolve() / 'bohra_defaults.json'}"
    with open(cfg_file, 'r') as f:
        cfg = json.load(f)
        return cfg["resource_levels"]
    
def _max_cpus(cpus:int) -> dict:

       
    avail = int(_get_cpu_limit_local(cpus=cpus))
    
    if cpus > avail:
        LOGGER.error(f"You requested {cpus} CPUs but only {avail} are available. Please your available resources and try again.")
        raise SystemExit
    
    return cpus

def _setup_config(default_config: dict, user_config: dict) -> dict:

    levs = True
    try:
        for level in default_config["levels"]:
            if level not in user_config["levels"]:
                levs = False
    except Exception as e:
        levs = False
        LOGGER.warning(f"Could not parse user config levels. Using default config levels: {', '.join(default_config['levels'])}. Error: {e}")
    if "profile" not in user_config:
        levs = False
        LOGGER.warning(f"Profile not found in user config. Using default profile: {default_config['profile']}.")
    if not levs:
        LOGGER.warning(f"One or more levels are missing from the user config. Using default config levels: {', '.join(default_config['levels'])}.")
        config = default_config
    else:
        config = user_config
    return config


def _check_config(user_config: str, default_config:dict) -> dict:

    if user_config != '':
        if _check_path(user_config):

            try:
                with open(user_config, 'r') as f:
                    user_cfg = json.load(f)
                    user_cfg = _setup_config(default_config=default_config, user_config=user_cfg)
                    return user_cfg
            except Exception as e:
                LOGGER.warning(f"Could not open profile config file {user_config}. Using the default profile 'lcl'. Error: {e}")
                return default_config
        else:
            LOGGER.warning(f"The profile config file {user_config} does not exist. Using the default profile 'lcl'.")
            return default_config
    return default_config

def _check_cfg_resources(cfg: dict, cpus:int) -> dict:

    for level in cfg["levels"]:
        if cpus < cfg["levels"][level]["cpus"]:
            cfg["levels"][level]["cpus"] = cpus
    return cfg

def _wrangle_config(config: dict) -> dict:

    res = {}
    res["profile_name"] = config["profile"]
    res["low"] = config["levels"]["low"]["cpus"]
    res["medium"] = config["levels"]["medium"]["cpus"]
    res["high"] = config["levels"]["high"]["cpus"]
    res["avail"] = os.cpu_count()
    return res

def _get_config(title: str, wd: str,user_config: str, cpus: int) -> str:

    _profile_cfg = _get_default_profile_cfg()

    LOGGER.info(f"Tyring to find your profile.")
    
    config = _check_config(user_config=user_config, default_config=_profile_cfg)
    config = _check_cfg_resources(cfg=config, cpus=cpus)
    tmpl = _wrangle_config(config=config)

    profile_config_path = _generate_nxf_config(title= title, res=tmpl, wd=wd)
    profile = config["profile"]

    return profile, profile_config_path