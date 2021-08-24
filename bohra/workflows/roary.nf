#!/usr/bin/env nextflow

include { ROARY } from './../modules/roary/main' 
include { ROARY2SVG } from './../modules/roary2svg/main' 

workflow RUN_ROARY {   

    take:
        gff
    main:
        ROARY ( gff )       
        ROARY2SVG ( ROARY.out.roary_csv )
    emit:
        svg = ROARY2SVG.out.pan_genome
        roary = ROARY.out.roary_summary

}

