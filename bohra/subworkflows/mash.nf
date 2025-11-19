#!/usr/bin/env nextflow

include { MASH_SKETCH } from './../modules/mash/sketch/main'
include { MASH_TRIANGLE } from './../modules/mash/triangle/main'
include { CSVTK_UNIQ } from './../modules/csvtk/main'


workflow RUN_MASH {

    take:
        sequences
    main:
        MASH_SKETCH ( sequences ) 
        sketches = MASH_SKETCH.out.sketch.map { cfg, sketch -> sketch }.collect()
        MASH_TRIANGLE ( sketches )
        versions = MASH_SKETCH.out.version.map { cfg, file -> file }.collect()
                                        .map { files -> tuple("version_mash", files) }
        CSVTK_UNIQ ( versions )
    emit:
        dists = MASH_TRIANGLE.out.dists
        version = CSVTK_UNIQ.out.collated
       

}

