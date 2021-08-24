[![CircleCI](https://circleci.com/gh/MDU-PHL/bohra.svg?style=svg&circle-token=530799cb0764519fc65966ab48bac7e0d02f3688)](https://circleci.com/gh/MDU-PHL/bohra)
[![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-370/)



# Bohra

Bohra is microbial genomics pipeline, designed predominantly for use in public health, but may also be useful in research settings. At a minimum the pipeline takes as input a tab-delimited file with the isolate IDs followed by the path to READ1 and READ2, a reference for alignment and a unique identifier, where reads are illumina reads (other platforms are not supported at this stage).

## Bohra has a new look! Welcome to Bohra-NF

**New features**

* Pipeline written in [Nextflow](https://www.nextflow.io)
* Default mode
    * Can be run on a single isolate (phylogenetic tree will not be generated if fewer than three sequences included in dataset)
    * MobSuite integration.
    * Updated [abriTAMR] with support for point mutations and virulence factors (beta).
* Roary with visualisation of pan-genome.
* Improved support for different computing environments.

**Accreditation**

* _snippy and snippy-core version 4.4.5 are NATA accredited for accurate detection of SNPs for reporting of genomic relationships at [MDU](https://biomedicalsciences.unimelb.edu.au/departments/microbiology-Immunology/research/services/microbiological-diagnostic-unit-public-health-laboratory#about-mdu-phl) Victoria Australia_ 
* _abritamr is accredited for detection of AMR genes at [MDU](https://biomedicalsciences.unimelb.edu.au/departments/microbiology-Immunology/research/services/microbiological-diagnostic-unit-public-health-laboratory#about-mdu-phl) Victoria Australia_ 
### Motivation

Bohra was inspired by Nullarbor (https://github.com/tseemann/nullarbor) to be used in public health microbiology labs for analysis of short reads from microbiological samples. The pipeline is written in [Nextflow](https://www.nextflow.io).

### Etymology

Bohra the name of an exinct species of tree kangaroo that lived on the Nullarbor plain in Australia. The name was chosen to reflect the fact that it will be predominantly used to build *trees*, relies on [*snippy*](https://github.com/tseemann/snippy) (named for a very famous kangaroo) and was inspired by [*nullarbor*](https://github.com/tseemann/nullarbor). 


## Pipeline

Bohra takes raw sequencing reads and produces a standalone html file for simple distribution of reports. An addtional file can be porvided with the paths to any assemblies that have already been generated. This is a helpful saver of time.

![Image](https://github.com/MDU-PHL/bohra/blob/master/workflow.png)

Bohra can be run in three modes
1. Preview
* Calculate mash-distances
* Build a mash-tree
* Report sequencing statistics
* Species identification (providing you have kraken2 database setup properly ;) )

2. Default SNPs, species ID, assemly, MLST, Resistome and annotation)
* Call variants
* Generate a phylogenetic tree
* Assemble 
* MLST
* Resistome
* Annotate
* Plasmid prediction
* Species identification

3. SNPs, Phylogeny, PanGenome and  Typing and Species Identification
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
* [Snippy (v4.4.5 recommended)](https://github.com/tseemann/snippy)
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
* [mob_suite](https://github.com/phac-nml/mob-suite)
* [csvtk](https://github.com/shenwei356/csvtk)


#### Check installation

Check that all dependencies are installed.

```
bohra --check
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
$ bohra -h

Bohra - a bacterial genomics pipeline - version 2.0.0

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  --check               Check that dependencies are installed correctly.
                        (default: False)
  --input_file INPUT_FILE, -i INPUT_FILE
                        Path to reads file, which is a tab-delimited with 3
                        columns <isolatename> <path_to_read1> <path_to_read2>.
                        REQUIRED (default: )
  --contigs CONTIGS, -c CONTIGS
                        Path to contigs file, which is a tab-delimited with 3
                        columns <isolatename> <path_to_contigs>. OPTIONAL if
                        you already have assemblies. (default: )
  --job_id JOB_ID, -j JOB_ID
                        Job ID, will be the name of the output directory
                        (default: )
  --reference REFERENCE, -r REFERENCE
                        Path to reference (.gbk or .fa) (default: )
  --mask MASK, -m MASK  Path to mask file if used (.bed) (default: )
  --abritamr_args {Acinetobacter_baumannii,Campylobacter,Enterococcus_faecalis,Enterococcus_faecium,Escherichia,Klebsiella,Salmonella,Staphylococcus_aureus,Staphylococcus_pseudintermedius,Streptococcus_agalactiae,Streptococcus_pneumoniae,Streptococcus_pyogenes,Vibrio_cholerae}
                        Set if you would like to use point mutations, please
                        provide a valid species. (default: )
  --kraken_db KRAKEN_DB, -k KRAKEN_DB
                        Path to DB for use with kraken2, if no DB present
                        speciation will not be performed. (default:
                        KRAKEN2_DEFAULT_DB)
  --pipeline {preview,default,all}, -p {preview,default,all}
                        The pipeline to run. `preview` - generates a rapid
                        tree using mash distances | `default` - runs snippy,
                        phylogenetic tree (if > 3 sequences), assemblies, mlst
                        and amr gene detection | `all` - same as default but
                        includes roary pangenome analysis (default: preview)
  --assembler {shovill,skesa,spades}, -a {shovill,skesa,spades}
                        Assembler to use. (default: spades)
  --cpus CPUS           Number of max CPU cores to run, will define how many
                        rules are run at a time (default: 16)
  --minaln MINALN, -ma MINALN
                        Minimum percent alignment. Isolates which do not align
                        to reference at this threshold will not be included in
                        core phylogeny. (default: 0)
  --minqual MINQUAL, -mq MINQUAL
                        Minimum Avg quality of reads (default: 0)
  --mincov MINCOV, -mc MINCOV
                        Minimum percent alignment. Isolates which do not have
                        average read coverage above this threshold will not be
                        included further analysis. (default: 0)
  --workdir WORKDIR, -w WORKDIR
                        The directory where Bohra will be run, default is
                        current directory.
  --force, -f           Add if you would like to force a complete restart of
                        the pipeline. All previous logs will be lost.
                        (default: False)
  --no_phylo            Set if you do NOT want to generate a phylogentic tree.
                        (default: False)
  --executor EXECUTOR   Type of queue
  --keep {Y,N}          If you are rerunning bohra over an exisiting directory
                        set --keep to 'Y' to archive report files - otherwise
                        previous reprot files will be removed. (default: N)

```


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



### Preview mode

`bohra` preview mode uses `mash` to calculate mash distances between isolates and generate a mash tree to rapidly identify outliers in your dataset or identify clades of interest for a more focused analysis.

```
bohra -i input.tab -r ref.fa -p preview -j job_id
```

### Default

Default mode will perform SNP detection, assemblies (if contigs file not provided), mlst, abritamr and phylogeny (unless `--no_phylo` or < 4 samples used as input).

```
bohra -i input.tab -c contigs.tab -r ref.fa -p default -j job_id
```


### All

In addition to the default mode above, `all` will also run `roary` . If less than 4 samples are provided, roary will not be run and pipeline will revert to `default`

```
bohra -i input.tab -c contigs.tab -r ref.fa -p all -j job_id
```

### Profile

You can provide a specific profile for running `bohra` on various computing platforms (see https://www.nextflow.io/docs/latest/executor.html for available platforms). You will need to specify a `--profile`, which MUST be the same as the name of the profile in your `--config` file. An example structure is provided in `example.config` 




