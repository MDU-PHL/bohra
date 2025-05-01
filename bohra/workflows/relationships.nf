#!/usr/bin/env nextflow

include { RUN_SNPS } from './../subworkflows/snippy' 
include { SNIPPY_CORE } from './../modules/snippy_core/main' 
include { SNP_DISTS } from './../modules/snp_dists/main' 
include { SNIPPY_QC } from './../modules/collation/main' 
include { GUBBINS } from './../modules/gubbins/main' 
include { SNIPPY_CLEAN } from './../modules/snippy_clean/main' 
include { IQTREE } from './../modules/iqtree/main' 
include { CSVTK_CONCAT } from './../modules/csvtk/main'


workflow RUN_SNPS {

    take:
        reads
        reference
    main:
        if ( params.modules.contains("snippy") ){
            SNIPPY ( reads.combine( reference ) )  
            SNIPPY_QC ( SNIPPY.out.aln )
            alns =  SNIPPY.out.aln.map { cfg, aln -> aln.getParent() }.collect()
        } else {
            alns = reads
        }
    // will import from subworkflows based on the modules selected
    // ska
    // snippy 
    // trees
    // dists
    // clusters
    // stats - where applicable

    
        
    // emit:
        
        // dists = SNP_DISTS.out.distances
        // clusters = SNP_DISTS.out.clusters
        // tree = tree
        // core_stats = all_core_stats
}

