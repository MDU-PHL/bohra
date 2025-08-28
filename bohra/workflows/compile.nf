#!/usr/bin/env nextflow

include { RUN_SMRYZR } from './../modules/datasmryzr/main'
include { CONCAT_FILES } from './../modules/utils/main'

workflow RUN_COMPILE {

    take:
        results
        versions

    main:
        versions = versions.collect()
                                .map { files -> tuple("versions", files) }
        CONCAT_FILES ( versions )

        results = CONCAT_FILES.out.collated.concat ( results ).collect()
        // println results.view()
        RUN_SMRYZR ( results )


    emit:
        report = RUN_SMRYZR.out.report

}