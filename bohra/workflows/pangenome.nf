#!/usr/bin/env nextflow

include { PANAROO } from './../modules/panaroo/main' 
include { ROARY2SVG } from './../modules/roary2svg/main' 

workflow RUN_PANAROO {   

    take:
        gff
    main:
        PANAROO ( gff )       
        ROARY2SVG ( PANAROO.out.pangenome_csv )
    emit:
        svg = ROARY2SVG.out.pan_genome
        roary = PANAROO.out.pangenome_summary

}

