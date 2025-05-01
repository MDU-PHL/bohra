#!/usr/bin/env nextflow

include { SHOVILL } from './../modules/shovill/main' 
include { SPADES } from './../modules/spades/main' params( options: [args2: "${params.spades_opt}"] )
include { SKESA } from './../modules/skesa/main' 
// include { COLLATE_ASM } from './../modules/collation/main'
// include { SEQKIT_STATS } from './../modules/seqkit/stats/main'
// include { PROKKA } from './../modules/prokka/main'
// include { CSVTK_CONCAT } from './../modules/csvtk/main'
workflow RUN_ASSEMBLE {   

    take:
        input
        
    main:
    reads = input.filter{ cfg,files -> cfg.modules.contains("assemble") && cfg.pe_reads == true }
    // contigs = input.filter{ cfg,files -> !cfg.modules.contains("assemble") && cfg.asm == true }

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
        } else if (params.assembler == 'flye'){
            println "Flye is not supported yet"
        } 
        contigs = contigs.map { cfg, files -> tuple([id: cfg.id, modules: cfg.modules, asm:true, pe_reads:false], files) }


        // SEQKIT_STATS ( contigs )
        // RUN_PROKKA ( contigs )
        // APS = RUN_PROKKA.out.prokka_txt.join( SEQKIT_STATS.out.stats )
        // COLLATE_ASM ( APS )
        // asm_stats = CSVTK_CONCAT ( COLLATE_ASM.out.assembly.map { cfg, asm -> asm }.collect().map { files -> tuple("assembly", files)} )
        
    emit:
        contigs
        // assembly_stats = asm_stats
        // gff = RUN_PROKKA.out.gff
        // prokka_txt = RUN_PROKKA.out.prokka_txt
        
}
