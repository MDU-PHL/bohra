# The home of bohra

![pretty](images/bohra_doc.png)

**Comprehensive sequence characterisation for microbial genomics**

## Introduction

`bohra` is microbial genomics pipeline, designed predominantly for use in public health, but may also be useful in research settings. It leverages existing high quality bioinformatics tools, to provide users with an easily accessible report of comprehensive analysis results of bacterial sequence data to for characterisation of single samples or for outbreak investigations or population studies. 

1. Quality assessment of the input data
2. Speciation and appropriate _in silico_ serotyping (where applicable).
3. MLST
4. Species relevant recovery of AMR mechanisms and inference of genomic AST/DST were available (_S. enterica_ and _M. tuberulosis_).
5. Plasmid information
6. Comparative analysis using a reference-free or reference-based appproaches.
7. Pangenome analysis.


The pipeline is designed to be flexible and modular, allowing for inputs from paired end fastq or assemblies, with direct support for ONT coming soon.

Stand alone html reports are generated for easy sharing and visualisation of the results.

## Installation

**Note the following instructions are for pre-release installation for `bohra` version 3**

`bohra` is a large pipeline with many dependencies, including databases required for speciation. Currently `bohra` is only available for installation with `conda`.

**1. Install conda (skip this step if you already have conda installed)**

If you do not already have `conda` installed, you can check out the documentation [here](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html). We recommend you install [`miniforge`](https://github.com/conda-forge/miniforge)


**2. Setup the `bohra` environment and install the pipeline co-ordinator**

1. Clone the `bohra` repository

```
git clone git@github.com:MDU-PHL/bohra.git
```


2. Create the environment 


```
cd bohra
conda env create -f environment.yml -n bohra
```

3. Install the `bohra` co-ordinator


```
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

It is highly recommended that you allow `bohra` to setup the required dependencies for the pipeline. 

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


## Etymology

The name 'bohra', is the name of an exinct species of tree kangaroo that lived on the Nullarbor plain in Australia was chosen to reflect the fact that it was originally developed to used to build trees, relies on [snippy](https://github.com/tseemann/snippy) (named for a very famous kangaroo) and was inspired by [nullarbor](https://github.com/tseemann/nullarbor).
