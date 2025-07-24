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

This step will setup conda environments in your `~/.conda` or `~/.miniconda3` (depending on how you have configured your `conda` installation). These environments will be prefixed 