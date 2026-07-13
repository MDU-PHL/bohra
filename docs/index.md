# The home of bohra

**Comprehensive sequence characterisation for microbial genomics - Documentation under contruction - thank you for your patience.**


The `bohra` pipeline is designed to be flexible and modular, allowing for inputs from paired-end fastq and/or assemblies, with direct support for ONT coming soon.

Stand alone html reports are generated for easy sharing and visualisation of the results.

Usage instructions can be found [here](usage/quickguide.md)

## Installation


### Recommended (conda or mamba)

**1. Install conda (skip this step if you already have conda installed)**

If you do not already have `conda` installed, you can check out the documentation [here](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html). We recommend you install [`miniforge`](https://github.com/conda-forge/miniforge)


**2. Create and install bohra**

```
mamba (or conda) create -n bohra -c bioconda bohra
```

**3. Install dependencies and setup databases**

```
conda activate bohra
bohra deps install
```

Details of database setup and other options for installation can be found [here](installation.md)

## Motivation


![pretty](images/bohra_doc.png)

`bohra` is microbial genomics pipeline, designed predominantly for use in public health, but may also be useful in research settings. It leverages existing high quality bioinformatics tools, to provide users with an easily accessible report of comprehensive analysis results of bacterial sequence data to for characterisation of single samples or for outbreak investigations or population studies. 

1. Quality assessment of the input data
2. Speciation and appropriate _in silico_ serotyping (where applicable).
3. MLST
4. Species relevant recovery of AMR mechanisms and inference of genomic AST/DST were available (_S. enterica_ and _M. tuberculosis_).
5. Plasmid information
6. Comparative analysis using a reference-free or reference-based appproaches.
7. Pangenome analysis.





## Etymology

The name 'bohra', is the name of an exinct species of tree kangaroo that lived on the Nullarbor plain in Australia and was chosen to reflect the fact that it was originally developed to used to build trees, relies on [snippy](https://github.com/tseemann/snippy) (named for a very famous kangaroo) and was inspired by [nullarbor](https://github.com/tseemann/nullarbor).
