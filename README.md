[![CI](https://github.com/MDU-PHL/bohra/actions/workflows/CI.yml/badge.svg)](https://github.com/MDU-PHL/bohra/actions/workflows/CI.yml)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/MDU-PHL/bohra)](https://github.com/MDU-PHL/bohra/releases/latest)
[![Conda Downloads](https://img.shields.io/conda/dn/bioconda/bohra)](https://anaconda.org/bioconda/bohra)
![Python](https://img.shields.io/badge/python-3.x-blue.svg)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# Bohra

<IMG SRC="docs/images/logo.svg" ALIGN="right" WIDTH="4ex" ALT="Bohra Logo">
Bohra is am extensive pipeline 
for taking genome sequences
(short reads or assemblies) 
and running common bioinformatics assays
across the, including
[genotyping, AMR detection, and phylogenetics](#workflow).

# Install
```
% conda create -n bohra -c bioconda bohra
% bohra deps install
% bohra --version
```

# Documentation

Read the [Bohra website](https://mdu-phl.github.io/bohra/)
to learn how to use all the availaile features.

# Workflow

<P><IMG SRC="workflow.png" ALIGN="left" WIDTH="75%" ALT="Bohra workflow"></P>

# Authors

* [Kristy Horna](https://github.com/kristyhoran)
* [Torsten Seemann](https://tseemann.github.io)
