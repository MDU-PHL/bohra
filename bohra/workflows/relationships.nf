#!/usr/bin/env nextflow

include { RUN_SNPS } from './../subworkflows/snippy' 
include { RUN_SKA } from './../subworkflows/ska'
include { MAKE_TREE } from './../subworkflows/trees'

workflow RELATIONSHIPS {

    take:
        sequences
        reference
    main:
        if ( params.modules.contains("snippy") ){
            RUN_SNPS ( sequences, reference )
            dists = RUN_SNPS.out.dists
            clusters = RUN_SNPS.out.clusters
            stats = RUN_SNPS.out.stats
            core_aln = RUN_SNPS.out.aln
            core_full_aln = RUN_SNPS.out.core_full_aln
        } 
        else if ( params.modules.contains("ska") ) {
            // add in a join to combine reads and asm
            RUN_SKA ( sequences )
            // dists = RUN_SKA.out.dists
            // clusters = RUN_SKA.out.clusters
            // stats = RUN_SKA.out.stats
            // core_aln = RUN_SKA.out.aln
            
        }
        if (params.tree_input == "snippy" & params.modules.contains("tree")) {
            
            MAKE_TREE ( RUN_SNPS.out.aln, RUN_SNPS.out.core_full_aln )
            tree = MAKE_TREE.out.newick
        } else if (params.tree_input == "ska" & params.modules.contains("tree")) {
            MAKE_TREE ( RUN_SKA.out.aln )
            tree = MAKE_TREE.out.newick
        } else {
            tree = Channel.empty()
        }
        
        
    // will import from subworkflows based on the modules selected
    // ska
    // snippy 
    // trees
    // dists
    // clusters
    // stats - where applicable

    
        
    // emit:
        
    //     dists 
    //     clusters
    //     tree 
    //     stats 
}

