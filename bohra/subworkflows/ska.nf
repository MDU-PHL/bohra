#!/usr/bin/env nextflow

include { SKA_ALIGN } from './../modules/ska/align/main' 
include { SKA_DISTANCE } from './../modules/ska/distance/main' 
include { SKA_MERGE } from './../modules/ska/merge/main' 
include {SKA_BUILD} from './../modules/ska/build/main'


workflow RUN_SKA {
    take:
        reads
        reference
    main:
        SKA_BUILD ( reads )
        SKA_MERGE ( SKA_BUILD.out.skf )
        merged_skf = SKA_MERGE.out.merged_skf
        SKA_ALIGN ( merged_skf )
        SKA_DISTANCE ( merged_skf )
        
    emit:
        dists = SKA_DISTANCE.out.distance_long
        matrix = SKA_DISTANCE.out.matrix
        aln = SKA_ALIGN.out.aln
}