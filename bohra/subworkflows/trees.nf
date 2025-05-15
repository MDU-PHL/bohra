#!/usr/bin/env nextflow

 
include { IQTREE } from './../modules/iqtree/main' 
include { VERYFASTTREE } from './../modules/veryfasttree/main' 
include { LOWER_TRIANGLE } from './../modules/lower_triangle/main'
include { QUICKTREE } from './../modules/quicktree/main' 


workflow MAKE_SNP_TREE {

    take:
        core_aln
        core_full_aln
        
    main:
        tree = Channel.empty()
        version = Channel.empty()
        if (params.modules.contains("iqtree") ){
            IQTREE ( core_aln, core_full_aln)
            tree = IQTREE.out.newick
            version = IQTREE.out.version
            // version = CSVTK_UNIQ ( versions )
        } else if ( params.modules.contains("veryfasttree") ){
            VERYFASTTREE ( core_aln )
            tree = VERYFASTTREE.out.newick
            version = VERYFASTTREE.out.version
        }
        tree = tree.ifEmpty( 'not_available' )
    emit:
        
        tree = tree
        version
}


workflow MAKE_DIST_TREE {

    take:
        dist
        
    main:
        tree = Channel.empty()

        if (params.modules.contains("ska") || params.modules.contains("snippy") ){
            LOWER_TRIANGLE ( dist )
            dist = LOWER_TRIANGLE.out.triangle

        }
        QUICKTREE ( dist )
        tree = QUICKTREE.out.newick

    
    emit:
        
        tree = tree
        version = QUICKTREE.out.version
}
