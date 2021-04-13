#!/usr/bin/env nextflow

include { SNIPPY } from './../modules/snippy/main' 
include { SNIPPY_QC } from './../modules/collation/main' 

workflow RUN_SNIPPY {   

    take:
        reads
    main:
        SNIPPY ( reads )       
        SNIPPY_QC ( SNIPPY.out.aln )
    emit:
        aln = SNIPPY.out.aln
        qual = SNIPPY_QC.out.snippy_qc

}

