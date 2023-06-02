#!/usr/bin/env nextflow

include { CSVTK_CONCAT } from './../modules/csvtk/main'
include { SHOVILL } from './../modules/shovill/main' 
include { SPADES } from './../modules/spades/main' addParams( options: [args2: "${params.spades_opt}"] )
include { SKESA } from './../modules/skesa/main' 
include { SEQKIT_STATS } from './../modules/seqkit/stats/main' addParams( options: [args2: 'contigs'] )
include { PROKKA } from './../modules/prokka/main'
include { COLLATE_ASM } from './../modules/collation/main'

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
        RUN_PROKKA ( contigs )
        APS = RUN_PROKKA.out.prokka_txt.join( SEQKIT_STATS.out.stats )
        COLLATE_ASM ( APS )
        asm_stats = CSVTK_CONCAT ( COLLATE_ASM.out.assembly.map { cfg, asm -> asm }.collect().map { files -> tuple("assembly", files)} )
        
    emit:
        contigs
        assembly_stats = asm_stats
        gff = RUN_PROKKA.out.gff
        prokka_txt = RUN_PROKKA.out.prokka_txt
        
}

workflow RUN_PROKKA {

    take:
        contigs
    main:
        PROKKA ( contigs )
    emit:
        gff = PROKKA.out.gff
        prokka_txt = PROKKA.out.prokka_txt

}

