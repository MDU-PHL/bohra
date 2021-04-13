#!/usr/bin/env nextflow

include { SHOVILL } from './../modules/shovill/main' addParams( options: [args2: params.assembler_threads] )
// include { SPADES } from './../modules/spades/main' addParams( options: [args2: params.assembler_threads] )
// include { SKESA } from './../modules/skesa/main' addParams( options: [args2: params.assembler_threads] )

workflow ASSEMBLE {   

    take:
        reads
    main:
    if ( params.assembler == 'shovill'){
        SHOVILL ( reads )   
        contigs = SHOVILL.out.contigs    
        } 
        // else if ( params.assembler == 'spades' ){
        // SPADES ( reads )
        // contigs = SPADES.out.contigs
        // } else if (params.assembler == 'skesa' ) {
        // SKESA ( reads )
        // contigs = SKESA.out.contigs
        // }

}