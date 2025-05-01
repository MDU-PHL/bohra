#!/usr/bin/env nextflow

 
include { IQTREE } from './../modules/iqtree/main' 
include { VERYFASTTREE } from './../modules/veryfasttree/main' 


workflow MAKE_TREE {

    take:
        core_aln
        core_full_aln
        
    main:
        
        if (params.modules.contains("iqtree") ){
            IQTREE ( core_aln, core_full_aln)
            tree = IQTREE.out.newick
        } else if ( params.modules.contains("veryfasttree") ){
            VERYFASTTREE ( core_aln )
            tree = VERYFASTTREE.out.newick
        }
        
    emit:
        
        tree = tree
}

