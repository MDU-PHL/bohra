#!/usr/bin/env nextflow
include { MOBSUITE } from './../modules/mobsuite/main'
include {CSVTK_CONCAT } from './../modules/csvtk/main'

workflow RUN_MOBSUITE {
    
    take:
        asm
        
        
    main:
        MOBSUITE ( asm )
        // println MOBSUITE.out.contig_report.view()
        // println MOBSUITE.out.mobs.view()

    emit:
        contig_report = MOBSUITE.out.contig_report
        mobs = MOBSUITE.out.mobs
        
}