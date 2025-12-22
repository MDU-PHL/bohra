
from bohra.launcher.Utils import CustomFormatter, _extract_tool_list
import logging
import subprocess
import pathlib
import os
import json

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


def _run_cmd(cmd:list)-> bool:

    """
    Run a command and log output.
    """
    LOGGER.info(f"Running command: {' '.join(cmd)}")
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding='utf-8')
    while process.poll() is None:
        l = process.stdout.readline().strip() # This blocks until it receives a newline.
        if len(l.split()) > 0:
            LOGGER.info(f"{l}")

    if process.returncode != 0:
        print(process.args)
        LOGGER.warning(f"Error running command: {' '.join(cmd)}")
        if len( process.stdout.read().strip()) > 0:
            LOGGER.warning(f"{process.stdout.read()}")
        return False
    else:
        LOGGER.info(f"Command completed successfully: {' '.join(cmd)}")
    return True


    
def _check_envs(cfg:dict)->bool:
    """
    Check if conda environments directory exists.
    """
    bohra_env_dir = f"{os.getenv('CONDA_PREFIX')}"
    target_envs_dir = f"{bohra_env_dir}/bohra_conda_envs"
    for env in cfg:

        for dep in cfg[env]:
            env_name = f"{target_envs_dir}/{env}"
            if not pathlib.Path(env_name).exists():
                LOGGER.warning(f"Conda environment {env} not found at {env_name}.")
                return False
            cmd = ["conda", "run", "-p", env_name]
            cmd.extend(dep.split())
            if not _run_cmd(cmd):
                LOGGER.critical(f"Dependency {env} not installed properly.")
                return False
    return True


def _install_envs(cfg:dict, envs_path:str, env:str="all",force_reinstall:bool=False)->bool:
    """
    Install conda environments.
    """
    # checking if mamba is installed
    installer = "conda"
    proc = subprocess.run(["which", "mamba"], capture_output=True, text=True)
    if proc.returncode == 0:
        installer = "mamba"
        LOGGER.info("Using mamba to install dependencies.")
    else:
        LOGGER.info("Using conda to install dependencies.")
    # check if envs path exists
    bohra_env_dir = f"{os.getenv('CONDA_PREFIX')}"
    if pathlib.Path(bohra_env_dir).exists():
        target_envs_dir = f"{bohra_env_dir}/bohra_conda_envs"
        LOGGER.info(f"Dependency environments will be installed in: {target_envs_dir}")
    else:
        LOGGER.critical("CONDA_PREFIX environment variable not set. Please activate a conda environment before installing dependencies.")
        raise SystemExit
    if env != "all":
        cfg = {env: cfg[env]}
    for env_name in cfg:
        if _check_envs(cfg={env_name: cfg[env_name]}) and not force_reinstall:
            LOGGER.info(f"Environment {env_name} already installed and force_reinstall is False. Skipping installation.")
            # continue
        else:
            LOGGER.info(f"Setting up environment: {env_name}")
            yml_file = f"{envs_path}/{env_name}.yml"
            if not pathlib.Path(yml_file).exists():
                LOGGER.critical(f"Environment file {yml_file} does not exist.")
                raise SystemExit
            cmd = [installer, "env", "create", "-f", yml_file, "-p", f"{target_envs_dir}/{env_name}"]
            if force_reinstall:
                cmd.append("--force")
            if not _run_cmd(cmd):
                return False
            else:
                if _check_envs( cfg={env_name: cfg[env_name]}):
                    LOGGER.info(f"Environment {env_name} installed and verified successfully.")
                else:
                    LOGGER.critical(f"Environment {env_name} failed verification after installation.")
                    return False
    LOGGER.info("All dependencies installed successfully and verified.")
    return True


def dependencies(_action:str = "install",
                       envs:str=f"{pathlib.Path(__file__).parent.parent}/environments",
                       config:str=f"{pathlib.Path(__file__).parent.parent}/config/dependencies.json",
                       kwargs: dict = {})->int:
    """
    Install bohra dependencies.
    """
    script_path = f"{pathlib.Path(__file__).parent}"
    dep_cfg = _extract_tool_list(config, tool=kwargs.get('tool', 'all'))
    actions = ['install', 'update', 'check']
    
    LOGGER.info(f"Will now try {_action} dependencies. Please be patient this may take some time!!... Maybe get coffee.")
    if _action not in actions:
        LOGGER.critical(f"Invalid action: {_action}. Must be one of {actions}.")
        raise SystemExit
    if _action == 'check':
        if not _check_envs(dep_cfg):
            LOGGER.critical("Some dependencies are missing or not installed properly.")
            return 1
        else:
            LOGGER.info("All dependencies are installed properly.")
            return 0
    elif _action in ['install', 'update']:
        force_reinstall = True if _action == 'update' else False
        env_to_install = kwargs.get('tool', 'all')
        if not _install_envs(dep_cfg, envs, env=env_to_install, force_reinstall=force_reinstall):
            LOGGER.critical("Error installing dependencies.")
            return 1
        else:
            return 0
            # else:
            #     LOGGER.critical("Some dependencies failed verification after installation.")
            #     return 1

            
        

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