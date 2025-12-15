# Installation

**Note the following instructions are for pre-release installation for `bohra` version 3**



**1. Install conda (skip this step if you already have conda installed)**

If you do not already have `conda` installed, you can check out the documentation [here](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html). We recommend you install [`miniforge`](https://github.com/conda-forge/miniforge)


**2. Setup the `bohra` environment and install the pipeline co-ordinator**

1. Clone the `bohra` repository

```
git clone git@github.com:MDU-PHL/bohra.git
```


2. Create the environment and install `bohra`


```
cd bohra
conda env create -f environment.yml -n bohra
pip3 install /path/to/the/repository_for/bohra
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
  generate-input  Generare input files for the Bohra pipeline.
  init-databases  Check that dependencies are installed correctly.
  install-deps    Install dependencies for Bohra - Highly recommended to...
  run             Run the Bohra pipeline.
  test            Check that bohra is installed correctly and runs as...
```

**3. Install dependencies and setup environment variables**


This step will setup conda environments under the `path/to/conda/envs/bohra`(depending on how you have configured your `conda` installation). 

Additionally, `bohra` depends upon either a correctly configure `kraken2` compatible database. The `bohra install-deps` command will optionally download the databases for your and also set the appropriate environment variables for you.

1. Activate your environment from step 2 above
```
conda activate bohra
```

2. Run the bohra dependency installation

```
bohra install-deps --setup_databases
```

The initial creation of the conda environments may take some time (~10 - 15 minutes). Once the environments are set up, `bohra` will try to set up your database environment variables. Although this is *not* essential it is **HIGHLY** recommended for ease of running and reproducibility. 

Please be aware of the following:

* If you elect to say no to setting up a `KRAKEN2_DEFAULT_DB` you will need to provide the `--kraken2_db` flag each time you run `bohra` if you wish to do speciation or use any features which depend on species (typing and AMR).

* If you are happy to use the `mlst` database that comes bundled with the `mlst` [tool](https://github.com/tseemann/mlst), then you can decline setting the `BOHRA_PUBMLST_DB` and `BOHRA_BLAST_DB`.


* If you already have a kraken2 database you can press 'n' when asked if you want to download. If you select 'y' please make sure that you have enough storage space for the databases. Some are quite large, up to 600 GB (half-terabyte).
