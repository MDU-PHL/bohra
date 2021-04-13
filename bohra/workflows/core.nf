#!/usr/bin/env nextflow

include { SNIPPY_CORE } from './../modules/snippy_core/main' 
// include { IQTREE } from './../modules/iqtree/main' 

workflow RUN_CORE {   

    take:
        alns
    main:
        SNIPPY_CORE ( alns )       
        // SNIPPY_QC ( SNIPPY.out.aln )
    // emit:
    //     aln = SNIPPY.out.aln
    //     qual = SNIPPY_QC.out.snippy_qc

}

