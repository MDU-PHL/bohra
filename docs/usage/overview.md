# Usage overview

Users have the option to supply paired-end fastq (support for ONT coming soon) and/or _de novo_ assemblies as inputs into the pipeline. 

- Paired-end fastq only - where you only have paired-end fastq `bohra` will generate _de novo_ assemblies if required.

- _de novo_ assemblies only - where you only have _de novo_ assemblies `bohra` will not use a reference-based approach for comparative analysis. However, reference-free comparative tools, `ska2` and `mash` are available.

- Both paired-end fastq and _de novo_ assemblies - in situations where you have _de novo_ assemblies already generate, you can supply both sequence types to `bohra`. It will use the supplied _de novo_ assemblies for any steps which require them, potentially saving you time.

Additionally, you can also supply the species value and any optional sample metadata that may be useful.

## Input file

`bohra` requires a single tab-delimited file as input, you can find examples [here] TO ADD.

| Column name | Description | Required? | 
| :---: | :---: | :---: | 
| Isolate | This is the name of the sequence or sample and will appear throughout `bohra` outputs. It must be unique.| Yes|
| r1 | The path to read 1 | If an assembly file is not supplied you must supply reads | 
| r2 | The path to read 2 | If an assembly file is not supplied you must supply reads |
| assembly | The path to the assembly for the isolate | If reads are not supplied you must supply an assembly file |
| Species_expected | The expected species of the sample or 'control'.  | No |

**Species column**

If you do not require speciation as part of the pipeline and already know the species, you can provide it here. Please note if no speciation is undertaken, `bohra` will use this value to undertake typing and AMR mechanisms/inferrence. If the species in this column is NOT accurate - unexpected results may occur.

**Controls**

If you include controls in an analysis you can add an `is_control` column to your input file and supply the type of control. Note that anything designated as a control will have single-sample analyses done (were possible) but will NOT be included in any comparative analyses.


**Tree Annotation and additional metadata**

You may also provide additional columns of relevant metadata in your input file. `bohra` will NOT do any data validation on these columns - that is up to the user. But any additional metadata provided will be visible on the tree provided in the report file and the summary table.

**`bohra` can generate the input file for you**

If you have 

* A table with a list of isolates and othe data (species or other metadata) (column 'Isolate' must be included) 

AND/OR

* Paths to your reads and/or contigs

`bohra` can generate the input file for you. 

```
bohra generate-input --isolate_ids <table_name>.txt --reads /path/to/reads --contigs /path/to/contigs --outname my_data.txt
```
This will generate a file which you can use as input into `bohra`.

**Note that on large file systems this may take a while**



## Pipelines

`bohra` is a flexible pipeline and allows users to customise the workflows used. Below is an overview of each workflow. More detail on tools and options for each workflow can be found [here](usage/running_bohra.md) and [here](usage/modules.md). Further explanations and detailed guides can be found [here](guides/overview.md)

| Pipelines | |
|:--- | :--- |
| **[basic](usage/basic)** | **[assembly](usage/assemble)**|
|**[amr and typing](usage/amr_typing)**| **[comparative analysis](usage/comparative)**|
|**[full](usage/full)** | **[tb](usage/tb)**|




## A note on databases

Many bioinformatics tools require the use of a database or data collection. Where possible and appropriate, `bohra` utilises the databases that come packaged with the tools being used in order to ensure expected behaviour and consistency. However, there are cases were the user will need to supply a database path. 

* `kraken2` - these databases are very large and not easily packaged for distribution with software. So the user will either need to have existing databases available or download them as part of the setup of `bohra`. This information is [here](../installation.md)

* `mlst` comes packaged with a collection of profiles and is ready to use. However, due to changes in licensing, the most up to date profiles cannot be included. As such if you have available to you a current mlst database that is configured for use with `mlst` you can set the `BOHRA_PUBMLST_DB` and `BOHRA_BLAST_DB` environment variables as part of the setup of `bohra`, as described [here](../installation.md)