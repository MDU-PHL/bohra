<!-- [![CircleCI](https://circleci.com/gh/MDU-PHL/bohra.svg?style=svg&circle-token=530799cb0764519fc65966ab48bac7e0d02f3688)](https://circleci.com/gh/MDU-PHL/bohra) -->
[![Python 3.7](https://img.shields.io/badge/python-3.9-blue.svg)](https://www.python.org/downloads/release/python-370/)


# Bohra

Bohra is microbial genomics pipeline, designed predominantly for use in public health, but may also be useful in research settings. At a minimum the pipeline takes as input a tab-delimited file with the isolate IDs followed by the path to READ1 and READ2, a reference for alignment and a unique identifier, where reads are illumina reads (other platforms are not supported at this stage).

* Pipeline written in [Nextflow](https://www.nextflow.io)
* Default mode
    * Can be run on a single isolate (phylogenetic tree will not be generated if fewer than three sequences included in dataset)
    * MobSuite integration.
    * Updated [abriTAMR](https://github.com/MDU-PHL/abritamr) with support for point mutations and virulence factors (beta).
* Panaroo with visualisation of pan-genome.
* Improved support for different computing environments.

## Recent changes to bohra
* Each bohra process is now run in its own conda environment. These will by default be included in your working directory.
* Alternatively you can install your own conda environments and provide the path to where these can be found (see help below).
* If you wish you can maintain the existing structure - where all the dependencies are installed in a single environment. BEWARE this can cause unexpected behaviour and even cause the pipeline to fail. If you wish to run bohra in this way please contact us for advice.

**Comming soon**

* Improved report structure
* Baby kraken
* Typing (_Salmonella_ spp., _Listeria monocytogenes_, _Neiserria_)
* Mtb AMR


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
* MLST
* Resistome
* Annotate
* Plasmid prediction
* Species identification
* Pan Genome

### Installation (with conda)

Bohra requires >=python3.9 and conda

See below for instructions on how to configure the databases for kraken2.

1. Set up conda - documentation for conda installation can be found [here](https://conda.io/en/latest/miniconda.html)
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

2. Create the conda environment
```
conda create -n <bohra_env_name> bohra python=3.9
conda activate <bohra_env_name>
```
You can also install `bohra` with 
```
pip3 install bohra
```

3. Ensure that you have a kraken2 database. Minikraken can obtained as follows
```
wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken2_v2_8GB_201904_UPDATE.tgz
tar -C $HOME -zxvf minikraken2_v2_8GB_201904_UPDATE.tgz
```
This will download and unzip the kraken2 DB. Other kraken2 DB are also available, you can find more information [here](https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads)

Once you have the DB downloaded you will have to create an environment variable called `KRAKEN2_DEFAULT_DB`. This can be done by adding the following to your `$HOME/.bashrc`
```
export KRAKEN_DEFAULT_DB=$HOME/minikraken2_v2_8GB_201904_UPDATE
```

#### Recommended conda environments (for the brave)

If you would like to use pre-installed conda envs - the following are suggested. You may also attempt to install all these dependencies into a single conda environment - however it is not recomended as there may be significant dependency clashes.

* [Snippy (v4.4.5 recommended)](https://github.com/tseemann/snippy)
* [Shovill (skesa and spades.py)](https://github.com/tseemann/shovill)
* [Panaroo](https://github.com/gtonkinhill/panaroo)
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



### Running bohra

**A word on resources**

* A minimum of 8 cpus is required for running `bohra`
* `bohra` has default settings for running at MDU (both on research and service systems). These ensure optimal running on these systems.
* Whether you are running `bohra` on MDU systems or not, you can override the default max cpus, by providing an upper limit for your job using the `--cpus` flag.
* If you are running on a queue or in the cloud you will need to provide an additional config file with details of your profile requirements.

#### Using CLI

**`bohra generate_input`**

`bohra` will generate an input file for you with the path to reads in the correct format, you will need to supply a path to where the reads can be found. You may also want to supply a file with a list of sample IDs (especially if there are more reads in the path provided than you wish to use). Please note `generate_input` assumes that the files are named as `<sample_ID_someother_stuff_{R1,R2}.f*q.gz>` and will output a file called `isolates.tab`

```
bohra generate_input -r path_to_reads
```

**`bohra run`**


**Input file**

The input file needs to be a tab-delimited file with three columns IsolateID, path to R1 and path to R2. 
```
Isolate-ID    /path/to/reads/R1.fq.gz    /path/to/reads/R2.fq.gz
```

In addition a file of contigs may also be provided if you have already generated assemblies you wish to use. This is also a tab-delimited file. 

```
Isolate-ID    /path/to/asm.fasta
```
**Reference**

The choice of reference is important for the accuracy of SNP detection and therefore the investigation of genomic relatedness. Appropriate references should be chosen following the guidelines below.
1. A closed reference from the same ST (where applicable) or a gold-standard reference (as may be used in M. tuberculosis).
2. A pacbio or nanopore assembly from MDU that is of the same type as the query dataset
3. A high quality de novo assembly of either an isolate in the dataset or an isolate of the same ST or type.

**Mask**

Phage masking is important for to prevent the inflation of SNPs that can be introduced by horizontal transfer as opposed to vertical transfer. For closed genomes or those that are publicly available `phaster-query.pl` can be used to identify regions for masking. If a denovo assembly is used the website `phaster.ca` can be used. Regions for masking should be provided in `.bed` format.



### Preview mode (default)

`bohra` preview mode uses `mash` to calculate mash distances between isolates and generate a mash tree to rapidly identify outliers in your dataset or identify clades of interest for a more focused analysis.

```
bohra run -i input.tab -r ref.fa -p preview -j job_id
```

### Default

Default mode will perform SNP detection, assemblies (if contigs file not provided), mlst, abritamr and phylogeny (unless `--no_phylo` or < 4 samples used as input).

```
bohra run -i input.tab -c contigs.tab -r ref.fa -p default -j job_id
```
In addition, if you would like `bohra` to output point mutations for AMR (based on `abritAMR` and AMRFinderplus), you can also add the `--abritamr_args` flag with one of the following species:
```
Neisseria
Acinetobacter_baumannii
Campylobacter
Enterococcus_faecalis
Enterococcus_faecium
Escherichia
Klebsiella
Salmonella
Staphylococcus_aureus
Staphylococcus_pseudintermedius
Streptococcus_agalactiae
Streptococcus_pneumoniae
Streptococcus_pyogenes,Vibrio_cholerae
```

### Plus pangenome

In addition to the default mode above, `all` will also run `panaroo` . If less than 4 samples are provided, `panaroo` will not be run and pipeline will revert to `default`

```
bohra run -i input.tab -c contigs.tab -r ref.fa -p pluspan -j job_id
```

### Profile

You can provide a specific profile for running `bohra` on various computing platforms (see https://www.nextflow.io/docs/latest/executor.html for available platforms). You will need to specify a `--profile`, which MUST be the same as the name of the profile in your `--config` file. An example structure is provided in `example.config` 




