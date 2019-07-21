---
title: 'Bohra: A pipeline for microbial genomics based on Snakemake'
tags:
  - Python
  - bioinformatics
  - microbial genomics
  - Snakemake
authors:
  - name: Kristy Horan
    orcid:
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: William Pitchers
    orcid: 0000-0000-0000-0000
    affiliation: 1
  - name: Mark Schultz
    orcid:
    affiliation: 1
  - name: Anders Goncalves da Silva
    orcid:
    affiliation: 1
  - name: Torsten Seemann
    orcid:
    affiliation: "1,2"
affiliations:
 - name: Microbiological Diagnostic Unit Public Health Laboratory, Department of Microbiology and Immunology, Peter Doherty Institute for Infection and Immunity, The University of Melbourne
   index: 1
 - name: Institution 2
   index: 2
date: 22 July 2019
bibliography: paper.bib
---

# Summary

Transition to genomics in public health needs well developed, maintainable
tools. As well as a commitment to maintenance in the long term. 

A crucial component of microbial genomics in public health is SNP detection
and phylogenetic analysis of bacterial whole-genome sequence data.

``Bohra`` is designed in Python and uses ``Snakemake`` workflow engine to 
drive the pipeline. 

It can run as single computer and on clusters.

It has singularity/docker support.

It is fast, analysing 900 TB isolates easily.

``Bohra`` returns an easily digestable and interactive HTML report.


# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this: ![Example figure.](figure.png)

# Acknowledgements

We acknowledge contributions ...

# References