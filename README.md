[![CircleCI](https://dl.circleci.com/status-badge/img/gh/MDU-PHL/bohra/tree/master.svg?style=svg)](https://dl.circleci.com/status-badge/redirect/gh/MDU-PHL/bohra/tree/master)
[![Python 3.9](https://img.shields.io/badge/python-3.9-blue.svg)](https://www.python.org/downloads/release/python-390/)


# Bohra

Bohra is microbial genomics pipeline, designed predominantly for use in public health, but may also be useful in research settings. At a minimum the pipeline takes as input a tab-delimited file with the isolate IDs followed by the path to READ1 and READ2, where reads are illumina reads (other platforms are not supported at this stage).

For detailed usage information please see our [wiki](https://github.com/MDU-PHL/bohra/wiki)

## Recent changes to bohra
* Install script to setup dependencies for you.
* babykraken download as part of the dependency installation.
* Addition of typers
    * [Kleborate](https://github.com/klebgenomics/Kleborate/wiki)
    * [stype](https://github.com/MDU-PHL/salmonella_typing) (NATA accredited ISO15189)
    * [meningotype](https://github.com/MDU-PHL/meningotype) (NATA accredited ISO15189)
    * [lissero](https://github.com/MDU-PHL/lissero) (NATA accredited ISO15189)
    * [ngmaster](https://github.com/MDU-PHL/ngmaster)
    * [ectyper](https://github.com/phac-nml/ecoli_serotyping)

**Comming soon**

* Improved report structure
* Mtb AMR


**Accreditation**

Many of the underlying tools of the bohra pipeline are NATA accredited by [MDU](https://biomedicalsciences.unimelb.edu.au/departments/microbiology-Immunology/research/services/) Victoria Australia (ISO1589).

* `snippy` and `snippy-core` version 4.4.5  
* `abritamr` 
* `stype`
* `meningotype`
* `lissero`

### Motivation

Bohra was inspired by Nullarbor (https://github.com/tseemann/nullarbor) to be used in public health microbiology labs for analysis of short reads from microbiological samples. The pipeline is written in [Nextflow](https://www.nextflow.io).

### Etymology

Bohra the name of an exinct species of tree kangaroo that lived on the Nullarbor plain in Australia. The name was chosen to reflect the fact that it will be predominantly used to build *trees*, relies on [*snippy*](https://github.com/tseemann/snippy) (named for a very famous kangaroo) and was inspired by [*nullarbor*](https://github.com/tseemann/nullarbor). 


