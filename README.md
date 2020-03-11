[![CircleCI](https://circleci.com/gh/MDU-PHL/bohra.svg?style=svg&circle-token=530799cb0764519fc65966ab48bac7e0d02f3688)](https://circleci.com/gh/MDU-PHL/bohra)
[![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-370/)

![Image](https://github.com/kristyhoran/bohra/blob/master/bohra.png)

# Bohra 

Bohra has a new look! 
* A new preview mode for 'sneak peak' at your dataset. 
* `rerun` has been deprecated. If you would like to rerun a job etc; use the run command.
* If rerunning a job, a new --keep flag to archive previous report files.
* Built in filtering features to automatically remove isolates with low average coverage or alignment to your reference.
* Provide a standardised config file for commonly used settings.

Bohra is microbial genomics pipeline, designed predominantly for use in public health, but may also be useful in research settings. The pipeline takes as input a tab-delimited file with the isolate IDs followed by the path to READ1 and READ2, a reference for alignment and a unique identifier, where reads are illumina reads (other platforms are not supported at this stage).

### Motivation

Bohra was inspired by Nullarbor (https://github.com/tseemann/nullarbor) to be used in public health microbiology labs for analysis of short reads from microbiological samples. The pipeline is written in [Snakemake](https://snakemake.readthedocs.io/en/stable/). 

### Etymology

Bohra the name of an exinct species of tree kangaroo that lived on the nullarbor. The name was chosen to reflect the fact that it will be predominantly used to build *trees*, relies on *snippy* (named for a very famous kangaroo) and was inspired by *nullarbor*. 


## Pipeline

Bohra takes raw sequencing reads and produces a standalone html file for simple distribution of reports.
![Image](https://github.com/MDU-PHL/bohra/blob/master/workflow.png)

Bohra can be run in four modes
1. Preview (DEFAULT)
* Calculate mash-distances
* Build a mash-tree

2. SNPs, species ID and Assembly based tools (MLST, Resistome and annotation)
* Clean reads
* Call variants
* Generate a phylogenetic tree
* Assemble
* Species identification
* MLST
* Resistome
* Annotate

4. SNPs, Phylogeny, PanGenome and  Typing and Species Identification
* Clean reads
* Call variants
* Generate a phylogenetic tree
* Assemble
* Species identification
* MLST
* Resistome
* Annotate
* Pan Genome

### Installation

Bohra requires >=python3.7

#### Conda (Highly recomended)

Installing bohra with conda will ensure that all dependencies are present. See below for instructions on how to configure the databases for kraken2.

Set up conda - documentation for conda installation can be found [here](https://conda.io/en/latest/miniconda.html)
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```
It is recomended that you set up a `bohra` environment
```
conda create -n <bohra_env_name> bohra
```
To use bohra
```
conda activate <bohra_env_name>
```
#### PyPi
If installing with `pip` you will need to ensure other dependencies are also installed.
```
pip3 install bohra
```
* [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
* [Snippy](https://github.com/tseemann/snippy)
* [Shovill (skesa and spades.py)](https://github.com/tseemann/shovill)
* [Roary](https://sanger-pathogens.github.io/Roary/)
* [Prokka](https://github.com/tseemann/prokka)
* [kraken2](https://ccb.jhu.edu/software/kraken/)
* [abritamr](https://github.com/MDU-PHL/abritamr)
* [mlst](https://github.com/tseemann/mlst)
* [iqtree](http://www.iqtree.org/)
* [seqtk](https://github.com/lh3/seqtk)
* [snp-dists](https://github.com/tseemann/snp-dists)
* [mash](https://github.com/lskatz/mashtree)


#### Check installation

Check that all dependencies are installed.

```
bohra check
```

*IMPORTANT*

In addition to installing kraken ensure that you have a kraken2 database. Minikraken can obtained as follows
```
wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken2_v2_8GB_201904_UPDATE.tgz
tar -C $HOME -zxvf minikraken2_v2_8GB_201904_UPDATE.tgz
```
This will download and unzip the kraken2 DB. Other kraken2 DB are also available, you can find more information [here](https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads)

Once you have the DB downloaded you will have to create an environment variable called `KRAKEN2_DEFAULT_DB`. This can be done by adding the following to your `$HOME/.bashrc`
```
export KRAKEN_DEFAULT_DB=$HOME/minikraken2_v2_8GB_201904_UPDATE
```

### Running bohra

#### Using CLI

```
bohra run -h
usage: bohra run [-h] [--input_file INPUT_FILE] [-S]
                 [--abritamr_singularity ABRITAMR_SINGULARITY]
                 [--job_id JOB_ID] [--reference REFERENCE] [--mask MASK]
                 [--kraken_db KRAKEN_DB] [--pipeline {preview,sa,all}]
                 [--assembler {shovill,skesa,spades}] [--cpus CPUS]
                 [--minaln MINALN] [--mincov MINCOV]
                 [--prefillpath PREFILLPATH] [-mdu] [-workdir WORKDIR]
                 [-resources RESOURCES] [-force] [-dry-run] [--cluster]
                 [--json JSON] [--queue QUEUE] [--gubbins] [--keep {Y,N}]

optional arguments:
  -h, --help            show this help message and exit
  --input_file INPUT_FILE, -i INPUT_FILE
                        Input file = tab-delimited with 3 columns
                        <isolatename> <path_to_read1> <path_to_read2>
                        (default: )
  -S, --use_singularity
                        Set if you would like to use singularity containers to
                        run bohra. (default: False)
  --abritamr_singularity ABRITAMR_SINGULARITY
                        The path to containers. If you want to use locally
                        stored contianers please pull from
                        dockerhub://mduphl/<toolname>. (default:
                        docker://mduphl/abritamr:v0.2.2)
  --job_id JOB_ID, -j JOB_ID
                        Job ID, will be the name of the output directory
                        (default: )
  --reference REFERENCE, -r REFERENCE
                        Path to reference (.gbk or .fa) (default: )
  --mask MASK, -m MASK  Path to mask file if used (.bed) (default: False)
  --kraken_db KRAKEN_DB, -k KRAKEN_DB
                        Path to DB for use with kraken2, if no DB present
                        speciation will not be performed. [env var:
                        KRAKEN2_DEFAULT_DB] (default: None)
  --pipeline {preview,sa,all}, -p {preview,sa,all}
                        The pipeline to run. Preview (--preview - default)
                        will calculate mash-distances and a mash-tree for
                        quick inspection of your dataset. SNPs and ASSEMBLIES
                        ('sa') will perform SNPs and ASSEMBLIES. ALL ('all')
                        will perform SNPS, ASSEMBLIES and ROARY for pan-genome
                        analysis (default: preview)
  --assembler {shovill,skesa,spades}, -a {shovill,skesa,spades}
                        Assembler to use. (default: shovill)
  --cpus CPUS, -c CPUS  Number of CPU cores to run, will define how many rules
                        are run at a time (default: 16)
  --minaln MINALN, -ma MINALN
                        Minimum percent alignment. Isolates which do not align
                        to reference at this threshold will not be included in
                        core phylogeny. (default: 80)
  --mincov MINCOV, -mc MINCOV
                        Minimum percent alignment. Isolates which do not have
                        average read coverage above this threshold will not be
                        included further analysis. (default: 40)
  --prefillpath PREFILLPATH, -pf PREFILLPATH
                        Path to existing assemblies - in the form
                        path_to_somewhere/isolatename/contigs.fa (default:
                        None)
  -mdu                  If running on MDU data (default: False)
  -workdir WORKDIR, -w WORKDIR
                        The directory where Bohra will be run, default is
                        current directory (default:
                        /home/khhor/dev/playground/bohra/20200218_/test_f)
  -resources RESOURCES, -s RESOURCES
                        Directory where templates are stored (default:
                        /home/khhor/dev/bohra/bohra/templates)
  -force, -f            Add if you would like to force a complete restart of
                        the pipeline. All previous logs will be lost.
                        (default: False)
  -dry-run, -n          If you would like to see a dry run of commands to be
                        executed. (default: False)
  --cluster             If you are running Bohra on a cluster. (default:
                        False)
  --json JSON           Path to cluster.json - required if --cluster is set
                        (default: )
  --queue QUEUE         Type of queue (sbatch or qsub currently supported) -
                        required if --cluster is set. (default: )
  --gubbins, -g         Set to use gubbins for recombination correction.
                        (default: False)
  --keep {Y,N}          If you are rerunning bohra over an exisiting directory
                        set --keep to 'Y' to archive report files - otherwise
                        previous reprot files will be removed. (default: N)

```

#### Using a config file 
Arguments that start with '--' (eg. --input_file) can also be set in a config file. Please put your config file in the working directory and name it 'bohra.conf'
Config file syntax allows: key=value, flag=true,
stuff=[a,b,c] (for details, see syntax at https://goo.gl/R74nmi). If an arg is
specified in more than one place, then commandline values override environment
variables which override config file values which override defaults.

### Preview mode

Bohra's preview mode (default mode) uses `mash` to calculate mash distances between isolates and generate a mash tree to rapidly identify outliers in your dataset or identify clades of interest for a more focused analysis.

### Set up

**Input file**

The input file needs to be a tab-delimited file with three columns IsolateID, path to R1 and path to R2. 
```
Isolate-ID    /path/to/reads/R1.fq.gz    /path/to/reads/R2.fq.gz
```

**Reference**

The choice of reference is important for the accuracy of SNP detection and therefore the investigation of genomic relatedness. Appropriate references should be chosen following the guidelines below.
1. A closed reference from the same ST (where applicable) or a gold-standard reference (as may be used in M. tuberculosis).
2. A pacbio or nanopore assembly from MDU that is of the same type as the query dataset
3. A high quality de novo assembly of either an isolate in the dataset or an isolate of the same ST or type.

**Mask**

Phage masking is important for to prevent the inflation of SNPs that can be introduced by horizontal transfer as opposed to vertical transfer. For closed genomes or those that are publicly available `phaster-query.pl` can be used to identify regions for masking. If a denovo assembly is used the website `phaster.ca` can be used. Regions for masking should be provided in `.bed` format.


### Run

**Minimal command to run in preview mode**

To use alternative modes, use `-p` with one of the following arguments.

`sa` phylogeny and assembly associated tools

`all` all functions (`sa` plus pan-genome analysis)


`bohra run -r path/to/reference -i path/to/inputfile -j unique_id -m path/to/maskfile (optional)`

### Running Bohra in a HPC environment
Bohra can be run in a HPC environment (currently only sbatch and qsub are supported). To do this some knowledge and experience in such environments is assumed. You will need to provide a file called `cluster.json`. This file will contain rule specifc and default settings for running the pipeline. An template is shown below (it is recommended that you use this template, settings have been established using a slurm queueing system), in addition you can see further documentation [here](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#cluster-configuration).

*example command*
```
bohra run -r <reference> -i <input_tab> -j <jobname> --cluster --json cluster.json --queue sbatch
```

```
bohra rerun --cluster --queue sbatch
```

*cluster.json*
```
{
    "__default__" :
    {
        "account" : "AccountName",
        "time" : "0-0:5:00",
        "cpus-per-task": "2",
        "partition" : "cloud",
        "mem" : "2G",
        "ntasks" : "4",
        "job" : "{rule}"
    },
    "snippy" :
    {
        
        "cpus-per-task" : "4",
        "time" : "0-0:5:00",
        "mem" : "8G" 
    },
    "assemble":
    {
        "cpus-per-task" : "4",
        "time" : "0-0:20:00",
        "mem" : "32G"

    },
    "kraken" :
    {
        "cpus-per-task" : "8",
        "time" : "0-0:20:00",
        "mem"  : "32G"
    },
    "run_iqtree_core" :
    {
        "time": "0-0:20:00"
    },
    "roary" : 
    {
        "cpus-per-task" : "36",
        "time" : "0-0:25:00",
        "mem": "8G"
    },
    "run_prokka" :
    {
        "cpus-per-task" : "8",
        "time" : "0-0:10:00"
    },
    "run_snippy_core" :
    {
        "time" : "0-0:15:00"
    }
}
```



