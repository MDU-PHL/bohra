#!/usr/bin/env nextflow

include { CSVTK_CONCAT } from './../modules/csvtk/main'
include { COLLATE_ASM } from './../modules/collation/main'
include { SHOVILL } from './../modules/shovill/main' addParams( options: [args2: params.assembler_threads] )
include { SPADES } from './../modules/spades/main' addParams( options: [args2: params.assembler_threads] )
include { SKESA } from './../modules/skesa/main' addParams( options: [args2: params.assembler_threads] )
include { SEQKIT_STATS } from './../modules/seqkit/stats/main' addParams( options: [args2: 'contigs'] )
include { ABRITAMR } from './../modules/resistome/main' addParams( options: [args2: 4] )
include { MLST } from './../modules/mlst/main' addParams( options: [args2: 4] )
include { PROKKA } from './../modules/prokka/main' addParams( options: [args2: params.prokka_threads] )


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
        PROKKA ( contigs )
        println SEQKIT_STATS.out.stats.view()
    emit:
        contigs
        assembly_stats = SEQKIT_STATS.out.stats
        resistome = ABRITAMR.out.matches
        mlst = MLST.out.mlst
        gff = PROKKA.out.gff
        prokka_txt = PROKKA.out.prokka_txt
}

workflow CONCAT_MLST {
    take:
        mlsts
    main:
        CSVTK_CONCAT ( mlsts )
    emit:
        collated_mlst = CSVTK_CONCAT.out.collated

}


workflow CONCAT_RESISTOMES {
    take:
        resistomes
    main:
        CSVTK_CONCAT ( resistomes )
    emit:
        collated_resistomes = CSVTK_CONCAT.out.collated

}


workflow CONCAT_ASM {
    take:
        asm
    main:
        CSVTK_CONCAT ( asm )
    emit:
        collated_assembly = CSVTK_CONCAT.out.collated

}


workflow COLLATE_ASM_PROKKA {
    take:
        asm
    main:
        COLLATE_ASM ( asm )
    emit:
        collated_asm = COLLATE_ASM.out.assembly

}