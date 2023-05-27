#!/usr/bin/env nextflow

include { ABRITAMR;COMBINE_AMR } from './../modules/resistome/main' addParams( options: [args2: 4] )
include { MLST } from './../modules/mlst/main' addParams( options: [args2: 4] )
// include { PROKKA } from './../modules/prokka/main' 

workflow TYPING {
    
    take:
        contigs
    main:
        ABRITAMR ( contigs )
        ab = ABRITAMR.out.abritamr_matches.join(ABRITAMR.out.abritamr_partials)
        COMBINE_AMR( ab )
        MLST ( contigs )
}