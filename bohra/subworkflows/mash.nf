#!/usr/bin/env nextflow

include { MASH_SKETCH } from './../modules/mash/sketch/main'
include { MASH_TRIANGLE } from './../modules/mash/triangle/main'



workflow RUN_MASH {

    take:
        sequences
    main:
        MASH_SKETCH ( sequences ) 
        sketches = MASH_SKETCH.out.sketch.map { cfg, sketch -> sketch }.collect()
        MASH_TRIANGLE ( sketches )
        
    emit:
        dists = MASH_TRIANGLE.out.dists
       

}

