#!/usr/bin/env nextflow

include { MASH_SKETCH } from './../modules/mash/sketch/main'
include { MASH_TRIANGLE } from './../modules/mash/triangle/main'
include { QUICKTREE } from './../modules/quicktree/main' 


workflow PREVIEW_NEWICK {

    take:
        preview
    main:
        MASH_SKETCH ( preview ) 
        sketches = MASH_SKETCH.out.sketch.map { cfg, sketch -> sketch }.collect()
        MASH_TRIANGLE ( sketches )
        QUICKTREE ( MASH_TRIANGLE.out.mash_distances )
    emit:
        nwk = QUICKTREE.out.preveiw_tree
       

}

