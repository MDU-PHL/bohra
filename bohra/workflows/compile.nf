#!/usr/bin/env nextflow

include { RUN_SMRYZR } from './../modules/datasmryzr/main'
include { CSVTK_UNIQ } from './../modules/csvtk/main'

workflow RUN_COMPILE {

    take:
        results
        versions

    main:
        versions = versions.collect()
                                .map { files -> tuple("versions", files) }
        CSVTK_UNIQ ( versions )

        results = CSVTK_UNIQ.out.collated.concat ( results ).collect()
        println results.view()
    //     RUN_SMRYZR ( results )


    // emit:
    //     report = RUN_SMRYZR.out.report

}