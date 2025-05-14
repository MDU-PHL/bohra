#!/usr/bin/env nextflow

include { RUN_SNPS } from './../subworkflows/snippy' 
include { RUN_SKA } from './../subworkflows/ska'
include { MAKE_TREE } from './../subworkflows/trees'

workflow TREE_GENERATION {

    take:
        dists
        core_aln
        core_full_aln
    main:
        tree = Channel.empty()
        if (params.tree_input == "alignment") {
            MAKE_SNP_TREE ( core_aln, core_full_aln.ifEmpty( 'no_full_aln' ) )
            tree = MAKE_TREE.out.newick

        } else if (params.tree_input == "distance" ) {
            MAKE_DIST_TREE ( dist )
            tree = MAKE_TREE.out.newick
        } 
        
        
    // will import from subworkflows based on the modules selected
    // ska
    // snippy 
    // trees
    // dists
    // clusters
    // stats - where applicable

    
        
    emit:
        
    //     dists 
    //     clusters
        tree 
    //     stats 
}

