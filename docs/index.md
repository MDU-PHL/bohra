# Welcome to datasmryzr documentation

![Python package](https://github.com/kristyhoran/datasmryzr/actions/workflows/python-package.yml/badge.svg)


Datasmryzr is a small command-line tool that is designed to collate and render pathogen genomics analysis results, such as trees, tables and matrices as a `.html` for sharing and interactive interrogation. 

## Installation

`datasmryzr` is a python package with simply dependencies. It is recommended to use a virtual or conda environment for installation.

```
<activate enviornment>
pip3 install git+https://github.com/kristyhoran/datasmryzr
```

## Usage

`datasmryzr` allows users to generate standalone `.html` files for visualisation, searching and sharing of genomic analysis results. It has been designed for pathogen genomics needs, but could also be used for other types of tabular and newick based data.

### Inputs

Any combination of tables, newick or vcf can be supplied for generation of the html report. However, there are some things to be aware of in order to get a good result.

* Tabular data (comma or tab-delimited files) can be used as input. 
    * All columns in input tables will be rendered - in the order that they are supplied.
    * File names will be used as menu labels in the `.html` file - so use sensible ones :wink:

* Pairwise matrix data (comma or tab-delimited files) can be supplied to generate distributions and heatmaps of distances.

* Newick file - a single newick file can be supplied per report.

* An annotation file can also be supplied if you want to be able to add colored annotations to the tree. This file can be one that has been already specified as tabular data or it can be an additional file.
    * The first column in the annotation file MUST contain names that match exactly to the tiplabels in the newick file.
    * If you only want a subset of values from the annotation file, you can supply these as a comma separated list in the command.
    * Only categorical variables can be visualised on the tree at the moment. If you want to display numerical values as categorical, you can provide this information in a configuration file (see below for details).

* If you would like to visualise the distribution of variants across the reference genome, you will need to also supply the vcf, with all samples in it (for example the core.vcf output from snippy), a reference genome and a mask file if one has been used.

#### Config

You can supply a configuration file, with comments and categorical columns (for tree annotation) using `--config` ([see below](#trees))

#### Utility options

As with most tools, you can supply various options to customise your report

* `--output` - this is the path to where you want your report saved. Defaults to current working directory.
* `--title` - you can give your report a title - this will appear at the top left-hand side of your document.
* `--description` - you can also add in a brief sentence to describe your report - this will appear as text below the title
* `--filename` - you can supply a custom name for the output file - defaults to `datasmryzr`
* `--background_color` - defaults to `#546d78` (blue grey). This is the color used in header, titles and bar graphs.
* `--font_color` - defaults to `#ffffff`.

## Cookbook

### Simple tables

You can generate a html with just tabular data. Any tablular data can be used, csv or tab-delimted files (support for .xlsx coming soon).

**Example command**

``` bash
datasmryzr --title 'A new report' --filename filename1.txt --filename filename2.csv --filename filename3.tsv
```

### Trees

Commonly pathogen genomics analyses involve the generation of a tree of some sort, which can be challenging to visualise and contextualise with other types of data. In order to generate a report with a tree and annotation, you will need a newick file and file with the data you wish to display on the tree. 

**Pro-tips**

* The first column of the annotation file must contain the tiplabels of the tree you wish to annotate. If a tiplabel is not present in the file, no annotation will be assigned to those tips.
* Only categorical values will be annotated on the tree. If you have numeric values that should be annotated as categorical on your tree you can supply a configuration file 

```json
{
    "comments":{
    "some_file_stub": "some comment"},
    "datatype": {
        "MLST":"input",
        "ST":"input"
    },
    "categorical_columns": [
        "MLST",
        "ST"
    ]

}
```

* Colours will be randomly assigned - we do not yet support custom colour schemes.

**Example command**

```bash
datasmryzr --tree tree.newick --annotate annotation_file.csv --annotate_cols cat1,cat4,cat5
```

### Core genome

Sometimes it is useful to see where SNPs are concentrated when you use core genome alignment to a reference genome. If you supply `datasmryzr` with 
1. vcf file, for example the `core.vcf` output from [snippy](https://www.google.com/search?client=firefox-b-d&q=snippy+github)
2. reference genome (fasta format or genbank) used for alignment
3. Alignment statisticts, for example the `core.txt` from [snippy](https://www.google.com/search?client=firefox-b-d&q=snippy+github).
4. A mask file in bed format (optional - advised if used in analysis as the masked regions will be colored in grey).

**Example command**

```bash
datasmryzr --filename core.txt --core-genome core.vcf -r ref.fa -m mask.bed
```

### Distance matrix

A very common question that is asked in microbial genomics is 'how far apart these things are', where things can be distances which represent SNPs, alleles or some other feature. If you supply `datasmryzr` with a distance matrix, you can generate heatmap and pairwise dsitributions plots in your report.

**Example command**

```bash
datasmryzr --distance-matrix distances.txt
```

## Exploring the html

### The Tree

NEED TO ADD IMAGE/VIDEO

By deafult if you have supplied a tree, this will be the first view that is loaded. 

* You can select which columns annotate onto the tree by selecting the `annotate` dropdown menu item and see the lgened by using the `toggle-legend` menu.
* On large trees you can select internal nodes to zoom to the that subtree.
* Trees can be exported as newick files or the current view exported as a png.

### Distances

NEED TO ADD IMAGE/VIDEO

If you have supplied a distance matrix, you will see a `Distance` menu, if you select this, you can choose a grap to display.
* Graph images can be saved by clicking on the three dots to the right of the image.
* Rolling over blocks on the heatmap will dispay tooltip with names of pairs and the distances.

### Core genome graph

The distribution of SNPs across the reference genome supplied is displayed in bins of SNPS per 5MB. You can use scroll to zoom in on a section of the genome and rollover the bars to see the number of SNPS at each position.

### Other tables

* Tablular data can be sorted and searched using the input boxes at the top of each column. 
* Adjust the width of a column by hovering over the edge of the column and dragging it left or right
* Download the tables for each table. As the creator - obviously this is not that useful to you - but can be useful for collaborators or others that you may share the data with.

