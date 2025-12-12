#!/usr/bin/env nextflow

include { RUN_SNPS } from './../subworkflows/snippy' 
include { RUN_SKA } from './../subworkflows/ska'
include { MAKE_SNP_TREE } from './../subworkflows/trees'
include { MAKE_DIST_TREE } from './../subworkflows/trees'
include { RUN_MASH } from './../subworkflows/mash'

workflow RELATIONSHIPS {

    take:
        sequences
        reference
    main:
        version = Channel.empty().ifEmpty( 'no_version' )
        tree_version = Channel.empty().ifEmpty( 'no_tree_version' )
        if ( params.modules.contains("snippy") ){
            RUN_SNPS ( sequences, reference )
            dists = RUN_SNPS.out.dists
            clusters = RUN_SNPS.out.clusters
            stats = RUN_SNPS.out.stats
            core_aln = RUN_SNPS.out.aln
            core_vcf = RUN_SNPS.out.core_vcf
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
            core_vcf = Channel.empty()
            version = RUN_SKA.out.versions
            
        } else if ( params.modules.contains("mash") ) {
            // add in a join to combine reads and asm
            println "Running mash"
            RUN_MASH ( sequences )
            dists = RUN_MASH.out.dists
            clusters = Channel.empty()
            stats = Channel.empty()
            core_aln = Channel.empty()
            core_full_aln = Channel.empty()
            core_vcf = Channel.empty()
            version = RUN_MASH.out.version
            
        }
        if (params.tree_input == "alignment" & params.modules.contains("tree")) {
            core_full_aln = core_full_aln.ifEmpty( 'no_full_aln' )
            MAKE_SNP_TREE ( core_aln, core_full_aln )
            tree = MAKE_SNP_TREE.out.tree
            tree_version = MAKE_SNP_TREE.out.version
        } else if (params.tree_input == "distance" & params.modules.contains("tree")) {
            MAKE_DIST_TREE ( dists )
            tree = MAKE_DIST_TREE.out.tree
            tree_version = MAKE_DIST_TREE.out.version
        } else {
            tree = Channel.empty()
        }
        
        
        
    emit:
        
        dists = dists.ifEmpty( 'no_results' )
        clusters = clusters.ifEmpty( 'no_results' )
        tree = tree.ifEmpty( 'no_results' )
        stats = stats.ifEmpty( 'no_results' )
        version = version.ifEmpty( 'no_version' )
        tree_version = tree_version.ifEmpty( 'no_tree_version' )
        core_vcf = core_vcf.ifEmpty( 'no_results' )
        
}