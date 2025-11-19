#!/usr/bin/env nextflow
include { MLST } from './../modules/mlst/main'
include {CSVTK_CONCAT;CSVTK_UNIQ } from './../modules/csvtk/main'
include { CONCAT_FILES } from './../modules/utils/main'
workflow RUN_MLST {
    
    take:
        asm
        
        
    main:
        MLST ( asm )
        mlst = MLST.out.mlst.map { cfg,file -> file}.collect().map { files -> tuple("mlst", files)}
        CONCAT_FILES ( mlst )
        collated_mlst = CONCAT_FILES.out.collated
        versions = MLST.out.version.map { cfg, file -> file }.collect()
                                        .map { files -> tuple("version_mlst", files) }
        CSVTK_UNIQ ( versions )
    emit:
        collated_mlst
        version = CSVTK_UNIQ.out.collated

}
