# Usage overview

## Input files

`bohra` requires a single tab-delimited file as input, you can find examples [here] TO ADD.

| Column name | Description | Required? | 
| :---: | :---: | :---: | 
| Isolate | This is the name of the sequence or sample and will appear throughout `bohra` outputs. It must be unique.| Yes|
| r1 | The path to read 1 | If an assembly file is not supplied you must supply reads | 
| r2 | The path to read 2 | If an assembly file is not supplied you must supply reads |
| assembly | The path to the assembly for the isolate | If reads are not supplied you must supply an assembly file |
| species | The expected species of the sample or 'control'.  | No |

### Species column

If you do not require speciation as part of the pipeline and already know the species, you can provide it here. Please note if no speciation is undertaken, `bohra` will use this value to undertake typing and AMR mechanisms/inferrence. If the species in this column is NOT accurate - unexpected results may occur. Furthermore if your analysis includes control sequences you can provide that information here (`control`) and the sequence will not be included in any comparative analysis.

### Annotation

Where you are undertaking a comparative analysis (`snippy`, `ska2`, `mash`) you may also provide additional columns of relevant metadata in your input file. `bohra` will do data validation on these columns - that is up to the user. But any additional metadata provided will be visible on the tree provided in the report file.

### `generate_input`

If you have 

* A table with a list of isolates and othe data (species or other metadata) (column 'Isolate' must be included) 

AND/OR

* Paths to your reads and/or contigs

`bohra` can generate the input file for you. 

```
bohra generate-input --isolate_ids <table_name>.tsv --reads /path/to/reads --contigs /path/to/contigs
```
This will generate a file called `bohra_input.tsv` which you can use as input into `bohra`.

**Note that on large file systems this may take a while**

## Speciation

By default `bohra` will use [`sylph`](https://sylph-docs.github.io) for speciation as it is quite quick and the database is relatively small, making it easier to install. However, this is ONLY available where your input are fastq files (as per `sylph` guidance). If you require speciation from assemblies to be undertaken, you will need to use `kraken2` as your speciation tool.

## A note on databases

Many bioinformatics tools require the use of a database or data collection. Where possible and appropriate, `bohra` utilises the databases that come packaged with the tools being used in order to ensure expected behaviour and consistency. However, there are cases were the user will need to supply a database path. 

* `kraken2` - these databases are very large and not easily packaged for distribution with software. So the user will either need to have existing databases available or download them as part of the setup of `bohra`. This information is [here](../installation.md)

* `sylph` - like the kraken2 databases this is a fairly large collection and needs to either already be available on the system or downloaded as part of the setup for `bohra`, as detailed [here](../installation.md)

* `mlst` comes packaged with a collection of profiles and is ready to use. However, due to changes in licensing, the most up to date profiles cannot be included. As such if you have available to you a current mlst database that is configured for use with `mlst` you can set the `BOHRA_PUBMLST_DB` and `BOHRA_BLAST_DB` environment variables as part of the setup of `bohra`, as described [here](../installation.md)


## Pipelines

### Basic

### Assemble

### Comparative

### AMR and typing

### TB

### Default

### Full

### _C. auris_ (COMING SOON)