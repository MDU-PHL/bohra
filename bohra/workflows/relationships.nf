#!/usr/bin/env nextflow

include { RUN_SNPS } from './../subworkflows/snippy' 
include { RUN_SKA } from './../subworkflows/ska'
include { MAKE_SNP_TREE } from './../subworkflows/trees'
include { MAKE_DIST_TREE } from './../subworkflows/trees'

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
            version = RUN_SNPS.out.version
        } 
        else if ( params.modules.contains("ska") ) {
            // add in a join to combine reads and asm
            RUN_SKA ( sequences )
            dists = RUN_SKA.out.dists
            clusters = RUN_SKA.out.clusters
            stats = RUN_SKA.out.stats
            core_aln = RUN_SKA.out.aln
            core_full_aln = Channel.empty()
            version = RUN_SKA.out.version
            
        }
        if (params.tree_input == "alignment" & params.modules.contains("tree")) {
            core_full_aln = core_full_aln.ifEmpty( 'no_full_aln' )
            MAKE_SNP_TREE ( core_aln, core_full_aln )
            tree = MAKE_TREE.out.newick
            tree_version = MAKE_TREE.out.version
        } else if (params.tree_input == "dist" & params.modules.contains("tree")) {
            MAKE_DIST_TREE ( dists )
            tree = MAKE_DIST_TREE.out.newick
            tree_version = MAKE_DIST_TREE.out.version
        } else {
            tree = Channel.empty()
        }
        
        
        
    emit:
        
        dists 
        clusters
        tree 
        stats 
        version
        tree_version
        
}