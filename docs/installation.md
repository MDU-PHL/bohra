# Installation

**Note the following instructions are for pre-release installation for `bohra` version 3**

`bohra` is a large pipeline with many dependencies, including databases required for speciation. Currently `bohra` is only available for installation with `conda`.

## Create the environment

### 1. Install conda

If you do not already have `conda` installed, you can check out the documentation [here](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html). We recommend you install `miniconda`

Below are instructions to install `miniconda` on a linux machine, however there are other distributions for Windows and MacOS [here](https://www.anaconda.com/docs/getting-started/miniconda/install)

1. Download and install miniconda 

```
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh
```
2. Restart your terminal and then 

```
source ~/miniconda3/bin/activate
```

3. Initialise `conda`

```
conda init --all
```

### 2. Setup the `bohra` environment and install the pipeline co-ordinator

1. You can download the `environment.yml` file [here](https://github.com/MDU-PHL/bohra/blob/rethink_structure/environment.yml). Please note that on the last line you may need to update the path to where your conda environments are stored. For example if you installed `miniconda` as above you should chnage this to line to:

```
prefix: ~/miniconda/envs/bohra-3pr
```
2. Create the environment 

```
conda env create -f environment.yml
```

3. Install the `bohra` co-ordinator

```
pip3 install git+https://github.com/MDU-PHL/bohra.git@rethink_structure
```
4. Test this has worked as expected
```
bohra --help
```
Should result in 

```
Usage: bohra [OPTIONS] COMMAND [ARGS]...

Options:
  --version  Show the version and exit.
  --help     Show this message and exit.

Commands:
  check           Check that dependencies are installed correctly.
  generate-input  Generare input files for the Bohra pipeline.
  install-deps    Install dependencies for Bohra - Highly recommended to...
  run             Run the Bohra pipeline.
  test            Check that bohra is installed correctly
```

### 3. Install dependencies and setup environment variables

It is highly recommended that you allow `bohra` to setup the required dependencies for the pipeline. 

This step will setup conda environments in your `~/.conda` or `~/.miniconda3` (depending on how you have configured your `conda` installation). These environments will be prefixed with the name of the environment that you have installed `bohra` into. For example if you used the `environment.yml` file in step 2 your prefix will be `bohra-3pr`. This will ensure consistency and prevent duplication of environments across a file system. 

This is also useful for public health users, were the versions of software and databases needs to be strictly controlled.

Additionally, `bohra` depends upon either a correctly configure `kraken2` compatible database OR a `sylph` compatible database. The `bohra install-deps` command will optionally download the databases for your and also set the appropriate environment variables for you.

1. Activate your environment from step 2 above
```
conda activate bohra-3pr
```

2. Run the bohra dependency installation

```
bohra install-deps
```

The initial creation of the conda environments may take some time. Once the environments are set up, `bohra` will try to set up your database environment variables. Although this is not essential it is HIGHLY recommended for ease of running and reproducibility. 

Please be aware of the following:

* If you elect to say no to setting up a `BOHRA_KRAKEN2_DB` or `BOHRA_SYLPH_DB`, you will need to provide EITHER the `--kraken2_db` OR a `--sylph_db` flag each time you run `bohra` if you wish to do speciation or use any features which depend on species (typing and AMR).

* If you are happy to use the `mlst` database that comes bundled with the `mlst` [tool](https://github.com/tseemann/mlst), then you can decline setting the `BOHRA_PUBMLST_DB` and `BOHRA_BLAST_DB`.

* If you are working in your own environments (`~/.conda` or `~/.miniconda`) you skip setting the `BOHRA_MOBSUITE_DB` environment variable. However, please note if you are setting up `bohra` in a share environment, you may run into permissions issues if `mobsuite` database requires updating at the time of running (this can happen) - so please ensure that you have a properly set up `mob_suite` database and provide this path when requested during set up.

* If you already have a kraken2 database and/or a sylph database you can press 'n' when asked if you want to download. If you select 'y' please make sure that you have enough storage space for the databases. Some are quite large, up to 600 GB (half-terabyte).

