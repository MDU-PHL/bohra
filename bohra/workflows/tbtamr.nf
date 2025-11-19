#!/usr/bin/env nextflow

include { TBTAMR } from './../modules/tbtamr/main'
include { CSVTK_CONCAT;CSVTK_UNIQ } from './../modules/csvtk/main'

workflow RUN_TBTAMR {

    take:
        reads

    main:
        TBTAMR ( reads )
        results = TBTAMR.out.tbtamr_txt.map { cfg, file -> file }.collect()
        results = results.map { files -> tuple("tbtamr", files) }
        CSVTK_CONCAT ( results )
        results = CSVTK_CONCAT.out.collated
        versions = TBTAMR.out.version.map { cfg, file -> file }.collect()
                                        .map { files -> tuple("version_tbtamr", files) }
        CSVTK_UNIQ ( versions )

    emit:
        results = results
        version = CSVTK_UNIQ.out.collated

}