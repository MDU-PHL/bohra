
from bohra.launcher.Utils import CustomFormatter, _extract_tool_list, _check_path
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
formatter = logging.Formatter('[%(levelname)s:%(asctime)s] %(message)s', datefmt='%Y-%m-%d %I:%M:%S %p') 
fh.setFormatter(formatter)
LOGGER.addHandler(ch) 
LOGGER.addHandler(fh)


def _run_cmd(cmd:list, check:bool=False)-> bool:

    """
    Run a command and log output.
    """
    if not check:
        LOGGER.info(f"Running command: {' '.join(cmd)}")
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding='utf-8')
    while process.poll() is None:
        l = process.stdout.readline().strip() # This blocks until it receives a newline.
        if len(l.split()) > 0 and not check:
            LOGGER.info(f"{l}")
    if process.returncode != 0:
        print(process.args)
        LOGGER.warning(f"Error running command: {' '.join(cmd)}")
        if len( process.stdout.read().strip()) > 0:
            LOGGER.warning(f"{process.stdout.read()}")
        return False
    else:
        if not check:
            LOGGER.info(f"Command completed successfully")
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
            if not _run_cmd(cmd, check=True):
                LOGGER.critical(f"Dependency {env} not installed properly.")
                return False
            else:
                LOGGER.info(f"{env} environment is found and {dep.split()[0]} appear to be installed properly.")
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
            tool = kwargs.get('tool', 'all')
            LOGGER.critical(f"There are missing or improperly installed dependencies. Please run 'bohra deps install' to install.")
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

def _get_db_config(config:str=f"{pathlib.Path(__file__).parent.parent}/config/databases.json")->dict:
    """
    Get database configuration from json file.
    """
    if not pathlib.Path(config).exists():
        LOGGER.critical(f"Database configuration file {config} does not exist.")
        raise SystemExit
    with open(config, 'r') as f:
        db_cfg = json.load(f)
    return db_cfg

def _construct_options_string(options:dict)->str:
    """
    Construct options string from dictionary.
    """
    options_list = []
    for val in options:
    
        line = f"{val['key']} - {val['name']} (size {val['size']})"
        options_list.append(line)
    return "\n".join(options_list)

def _check_databases(db_install:bool=False,
                    database_path:str="",
                    config:str=f"{pathlib.Path(__file__).parent.parent}/config/databases.json")->int:
    """
    Check that required databases are installed.
    """
    script_path = f"{pathlib.Path(__file__).parent}"
    LOGGER.info(f"Will now check databases required for Bohra. Please be patient this may take some time!!... Maybe get coffee.")
    db_check_cmd = 'get' if db_install else 'check'
    db_cfg = _get_db_config(config=config)
    
    for essential in db_cfg["essential"]:
        var = os.getenv(essential)
        
        if var is None or not pathlib.Path(var).exists():
            LOGGER.warning(f"Essential database {essential} not found or environment variable not set.")
            if var == "KRAKEN2_DB_PATH":
                LOGGER.info(f"You will need to provide the path to the Kraken2 database for Bohra run speciation and other species-specific analyses.")
            if not db_install:
                LOGGER.warning(f"If you wish to install missing databases, please run 'bohra init_databases --setup_databases --database_path <path>'.")
            else:
                LOGGER.info(f"You will need to select one of the following options to setup the database:")
                options_str = _construct_options_string(db_cfg["essential"][essential]["urls"])
                db_choice = input(f"Select an option to setup {essential}:\n{options_str}\nEnter option key: ")
                selected_url = None
                for option in db_cfg["essential"][essential]["urls"]:
                    if str(option["key"]) == db_choice:
                        selected_url = option["url"]
                        target = option["target"]
                        break
                if selected_url is None and db_choice != "0":
                    LOGGER.critical(f"Invalid option selected: {db_choice}. Exiting.")
                    raise SystemExit
                elif db_choice == "0":
                    existing_path = input(f"Please provide the existing path to the {essential} database: ")
                    if _check_path(existing_path):
                        new_essential_path = existing_path
                        # LOGGER.info(f"Environment variable {essential} set to {existing_path}.")
                        continue
                    else:
                        LOGGER.critical(f"Provided path {existing_path} is invalid. Exiting.")
                        raise SystemExit
                if database_path == "":
                    database_path = input(f"Please provide the path to download and setup the {essential} database: ")
                if _check_path(database_path):
                    cmd = f"wget --continue -O {database_path}/{target} {selected_url}"
                    LOGGER.info(f"Running : {cmd}")
                    if not _run_cmd(cmd.split()):
                        LOGGER.critical(f"Error downloading database from {selected_url}. Exiting.")
                        raise SystemExit
                    else:
                        new_essential_path = f"{database_path}/{pathlib.Path(target).name.split('.')[0]}"
                        LOGGER.info(f"Uncompressing database {target}. This may take some time...")
                        cmd = f"tar -xvzf {database_path}/{target} -C {database_path}"
                        LOGGER.info(f"Running : {cmd}")
                        if not _run_cmd(cmd.split()):
                            LOGGER.critical(f"Error uncompressing database {target}. Exiting.")
                            raise SystemExit
                        LOGGER.info(f"Database {essential} downloaded.")
                LOGGER.info(f"Setting environment variable {essential} to {new_essential_path}.")
                cmd = f"conda env config vars set {essential}={new_essential_path} && conda deactivate && conda activate {os.getenv('CONDA_PREFIX')}"
                LOGGER.info(f"Running : {cmd}")
                # if not _run_cmd(cmd.split()):
                #     LOGGER.critical(f"Error setting environment variable {essential}. Exiting.")
                #     raise SystemExit
        else:
            LOGGER.info(f"Database {essential} found at {var}.")
    for other in db_cfg["non_essential"]:
        var = os.getenv(other)
        if var is None or not pathlib.Path(var).exists():
            LOGGER.warning(f"Optional database {other} not found or environment variable not set.")
            LOGGER.warning(f"Bohra will use the default databases that come with the tools for analyses requiring {other}. Otherwise you may be able to provide this at runtime.")
        else:
            LOGGER.info(f"Database {other} found at {var}.")
   