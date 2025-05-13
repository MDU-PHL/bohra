#!/usr/bin/env nextflow
include { MLST } from './../modules/mlst/main'
include {CSVTK_CONCAT } from './../modules/csvtk/main'

workflow RUN_MLST {
    
    take:
        asm
        
        
    main:
        MLST ( asm )
        mlst = MLST.out.mlst.map { cfg,file -> file}.collect().map { files -> tuple("mlst", files)}
        CSVTK_CONCAT ( mlst )
        collated_mlst = CSVTK_CONCAT.out.collated
    emit:
        collated_mlst

}
