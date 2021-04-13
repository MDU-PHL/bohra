#!/usr/bin/env nextflow

include { SHOVILL } from './../modules/shovill/main' addParams( options: [args2: params.assembler_threads] )
include { SPADES } from './../modules/spades/main' addParams( options: [args2: params.assembler_threads] )
include { SKESA } from './../modules/skesa/main' addParams( options: [args2: params.assembler_threads] )
include { SEQKIT_STATS } from './../modules/seqkit/stats/main' addParams( options: [args2: 'contigs'] )
include { ABRITAMR } from './../modules/resistome/main' addParams( options: [args2: 4] )
include { MLST } from './../modules/mlst/main' addParams( options: [args2: 4] )

workflow RUN_ASSEMBLE {   

    take:
        reads
    main:
    if ( params.assembler == 'shovill'){
        SHOVILL ( reads )   
        contigs = SHOVILL.out.contigs    
        } 
        else if ( params.assembler == 'spades' ){
        SPADES ( reads )
        contigs = SPADES.out.contigs
        } else if (params.assembler == 'skesa' ) {
        SKESA ( reads )
        contigs = SKESA.out.contigs
        }
        SEQKIT_STATS ( contigs )
        ABRITAMR ( contigs )
        MLST ( contigs )
    
    emit:
        contigs
        assembly_stats = SEQKIT_STATS.out.stats
        resistome = ABRITAMR.out.matches
        mlst = MLST.out.mlst
}