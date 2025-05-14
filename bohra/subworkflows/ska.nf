#!/usr/bin/env nextflow

include { SKA_ALIGN } from './../modules/ska/align/main' 
include { SKA_DISTANCE } from './../modules/ska/distance/main' 
include { SKA_MERGE } from './../modules/ska/merge/main' 
include {SKA_BUILD} from './../modules/ska/skf/main'
include { SNP_CLUSTER } from './../modules/cluster/main'

workflow RUN_SKA {
    take:
        sequence

        
    main:
        // println sequence.view()
        SKA_BUILD ( sequence )

        skf =  SKA_BUILD.out.skf.map { cfg, skf -> skf }.collect()
        // println skf.view()
        SKA_MERGE ( skf )
        merged_skf = SKA_MERGE.out.merged_skf
        println merged_skf.view()
        SKA_ALIGN ( merged_skf )
        SKA_DISTANCE ( merged_skf )
        SNP_CLUSTER ( SKA_DISTANCE.out.matrix )
    emit:
        stats = SKA_DISTANCE.out.distance_long
        dists = SKA_DISTANCE.out.matrix
        aln = SKA_ALIGN.out.aln
        clusters = SNP_CLUSTER.out.clusters
}