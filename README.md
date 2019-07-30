[![CircleCI](https://circleci.com/gh/MDU-PHL/bohra.svg?style=svg&circle-token=530799cb0764519fc65966ab48bac7e0d02f3688)](https://circleci.com/gh/MDU-PHL/bohra)
[![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-370/)

![Image](https://github.com/kristyhoran/bohra/blob/master/bohra.png)

# Bohra 

Bohra is microbial genomics pipeline, designed predominantly for use in public health, but may also be useful in research settings. The pipeline takes as input a tab-delimited file with the isolate IDs followed by the path to READ1 and READ2, a reference for alignment and a unique identifier, where reads are illumina paired end reads (other platforms are not supported).

### Motivation

Bohra was inspired by Nullarbor (https://github.com/tseemann/nullarbor) to be used in public health microbiology labs for analysis of short reads from microbiological samples. The pipeline is written in [Snakemake](https://snakemake.readthedocs.io/en/stable/). 

### Etymology

Bohra the name of an exinct species of tree kangaroo that lived on the nullarbor. The name was chosen to reflect the fact that it will be predominantly used to build *trees*, relies on *snippy* (named for a very famous kangaroo) and was inspired by *nullarbor*. 


## Pipeline

Bohra takes raw sequencing reads and produces a standalone html file for simple distribution of reports.
![Image](https://github.com/MDU-PHL/bohra/blob/master/workflow.png)

Bohra can be run in three modes
1. SNPs and Phylogeny
* Clean reads
* Call variants
* Generate a phylogenetic tree

2. SNPs, Phylogeny, Typing, Annotation and Species Identification (DEFAULT)
* Clean reads
* Call variants
* Generate a phylogenetic tree
* Assemble
* Species identification
* MLST
* Resistome
* Annotate

3. SNPs, Phylogeny, PanGenome and  Typing and Species Identification
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

#### Dependencies

Bohra requires >=python3.6 
```
pip3 install bohra
```
A conda recipe for bohra will follow shortly! But for now you will need to have the following dependencies installed on your system

* [Snippy](https://github.com/tseemann/snippy)
* [Shovill (skesa and spades.py)](https://github.com/tseemann/shovill)
* [Roary](https://sanger-pathogens.github.io/Roary/)
* [Prokka](https://github.com/tseemann/prokka)
* [kraken2](https://ccb.jhu.edu/software/kraken/)
* [abricate](https://github.com/tseemann/abricate)
* [mlst](https://github.com/tseemann/mlst)
* [iqtree](http://www.iqtree.org/)
* [seqtk](https://github.com/lh3/seqtk)
* [snp-dists](https://github.com/tseemann/snp-dists)

Bohra can be run in two modes `run` for an initial analysis and `rerun` for a re-analysis. A `.html` report is generated allowing for the visualisation of tree and examination of the dataset to provide insights that may be useful in interpretation of the results.

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

**Minimal command**

`bohra run -r path/to/reference -i path/to/inputfile -j unique_id -m path/to/maskfile (optional)`

```bohra run -h
usage: bohra run [-h] [--input_file INPUT_FILE] [--job_id JOB_ID]
                 [--reference REFERENCE] [--mask MASK]
                 [--pipeline {sa,s,a,all}]
                 [--assembler {shovill,skesa,spades}] [--cpus CPUS]
                 [--minaln MINALN] [--prefillpath PREFILLPATH] [--mdu MDU]
                 [--workdir WORKDIR] [--resources RESOURCES] [--force]
                 [--dryrun] [--gubbins]

optional arguments:
  -h, --help            show this help message and exit
  --input_file INPUT_FILE, -i INPUT_FILE
                        Input file = tab-delimited with 3 columns
                        <isolatename> <path_to_read1> <path_to_read2>
                        (default: )
  --job_id JOB_ID, -j JOB_ID
                        Job ID, will be the name of the output directory
                        (default: )
  --reference REFERENCE, -r REFERENCE
                        Path to reference (.gbk or .fa) (default: )
  --mask MASK, -m MASK  Path to mask file if used (.bed) (default: False)
  --pipeline {sa,s,a,all}, -p {sa,s,a,all}
                        The pipeline to run. SNPS ('s') will call SNPs and
                        generate phylogeny, ASSEMBLIES ('a') will generate
                        assemblies and perform mlst and species identification
                        using kraken2, SNPs and ASSEMBLIES ('sa' - default)
                        will perform SNPs and ASSEMBLIES. ALL ('all') will
                        perform SNPS, ASSEMBLIES and ROARY for pan-genome
                        analysis (default: sa)
  --assembler {shovill,skesa,spades}, -a {shovill,skesa,spades}
                        Assembler to use. (default: shovill)
  --cpus CPUS, -c CPUS  Number of CPU cores to run, will define how many rules
                        are run at a time (default: 36)
  --minaln MINALN, -ma MINALN
                        Minimum percent alignment (default: 0)
  --prefillpath PREFILLPATH, -pf PREFILLPATH
                        Path to existing assemblies - in the form
                        path_to_somewhere/isolatename/contigs.fa (default:
                        None)
  --mdu MDU             If running on MDU data (default: True)
  --workdir WORKDIR, -w WORKDIR
                        Working directory, default is current directory
                        (default: /home/khhor)
  --resources RESOURCES, -s RESOURCES
                        Directory where templates are stored (default:
                        /home/khhor/dev/bohra/bohra/templates)
  --force, -f           Add if you would like to force a complete restart of
                        the pipeline. All previous logs will be lost.
                        (default: False)
  --dryrun, -n          If you would like to see a dry run of commands to be
                        executed. (default: False)
  --gubbins, -g         If you would like to run gubbins. NOT IN USE YET -
                        PLEASE DO NOT USE (default: False)
```

### Rerun

A rerun may be performed if changes to the reference and/or mask file are needed. In addition, if isolates need to be removed or added to the analysis. 
The following behaviour on a rerun should be expected;
* New reference will result in calling of snps in all isolates of the analysis
* If the reference is unchanged SNPs will only be called on new isolates
* Determination of core alignment, distances and generation of trees will occur for every rerun

`-r` and `-m` are only required if these are to be different to the previous run. If not Bohra will detect and use the previous reference and mask files. Also changes to the isolates included should be made to the input file used in the original run. New isolates can be added to the bottom of the input file and prefixing an isolate with `#` will remove it from the analysis.

**Minimal command**

`bohra rerun`

```usage: bohra rerun [-h] [--reference REFERENCE] [--mask MASK] [--cpus CPUS]
                   [--workdir WORKDIR] [--resources RESOURCES] [--dryrun]
                   [--gubbins] [--keep]

optional arguments:
  -h, --help            show this help message and exit
  --reference REFERENCE, -r REFERENCE
                        Path to reference (.gbk or .fa) (default: )
  --mask MASK, -m MASK  Path to mask file if used (.bed) (default: )
  --cpus CPUS, -c CPUS  Number of CPU cores to run, will define how many rules
                        are run at a time (default: 36)
  --workdir WORKDIR, -w WORKDIR
                        Working directory, default is current directory
                        (default: /home/khhor)
  --resources RESOURCES, -s RESOURCES
                        Directory where templates are stored (default:
                        /home/khhor/dev/bohra/bohra/templates)
  --dryrun, -n          If you would like to see a dry run of commands to be
                        executed. (default: False)
  --gubbins, -g         If you would like to run gubbins. NOT IN USE YET -
                        PLEASE DO NOT USE (default: False)
  --keep, -k            Keep report from previous run (default: False)```

