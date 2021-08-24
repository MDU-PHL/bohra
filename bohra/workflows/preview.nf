#!/usr/bin/env nextflow

include { MASH_TRIANGLE } from './../modules/mash/triangle/main'
include { QUICKTREE } from './../modules/quicktree/main' 


workflow PREVIEW_NEWICK {

    take:
        sketches
        // seqkit
        // seqtk_stats
    main:
         MASH_TRIANGLE ( sketches )
         QUICKTREE ( MASH_TRIANGLE.out.mash_distances )
    emit:
        nwk = QUICKTREE.out.preveiw_tree

}

