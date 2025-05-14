#!/usr/bin/env nextflow

include { RUN_SNPS } from './../subworkflows/snippy' 
include { RUN_SKA } from './../subworkflows/ska'
include { RUN_MASH } from './../subworkflows/mash'
include { MAKE_SNP_TREE;MAKE_DIST_TREE } from './../subworkflows/trees'
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
            dists = RUN_SKA.out.dists
            clusters = RUN_SKA.out.clusters
            stats = RUN_SKA.out.stats
            core_aln = RUN_SKA.out.aln
            core_full_aln = Channel.empty().ifEmpty( 'no_full_aln' )
            
        } else if ( params.modules.contains("mash")){
            RUN_MASH ( sequences )
            dists = RUN_MASH.out.dists
            clusters = Channel.empty().ifEmpty( 'no_clusters' )
            stats = Channel.empty().ifEmpty( 'no_stats' )
            core_aln = Channel.empty().ifEmpty( 'no_core_aln' )
            core_full_aln = Channel.empty().ifEmpty( 'no_full_aln' )
        }
        
        tree = Channel.empty()
        if (params.modules.contains("tree")){
            if (params.tree_input == "alignment") {
                MAKE_SNP_TREE ( core_aln, core_full_aln )
                tree = MAKE_SNP_TREE.out.tree

            } else if (params.tree_input == "distance" ) {
                MAKE_DIST_TREE ( dists )
                tree = MAKE_DIST_TREE.out.tree
            }
        } else {
            tree = Channel.empty().ifEmpty( 'no_tree' )
        }
    
        
    emit:
        core_full_aln
        core_aln
        dists 
        clusters
        stats 
        tree
}

