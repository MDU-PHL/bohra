#!/usr/bin/env nextflow
include { MOBSUITE } from './../modules/mobsuite/main'
include {CSVTK_CONCAT;CSVTK_UNIQ } from './../modules/csvtk/main'

workflow RUN_MOBSUITE {
    
    take:
        asm
        
        
    main:
        MOBSUITE ( asm )
        // println MOBSUITE.out.contig_report.view()
        // println MOBSUITE.out.mobs.view()
        versions = MOBSUITE.out.version.map { cfg, file -> file }.collect()
                                        .map { files -> tuple("version_plasmid", files) }
        CSVTK_UNIQ ( versions )
    emit:
        contig_report = MOBSUITE.out.contig_report
        mobs = MOBSUITE.out.mobs
        version = CSVTK_UNIQ.out.collated
        
}