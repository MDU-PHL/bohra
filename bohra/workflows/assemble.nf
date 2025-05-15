#!/usr/bin/env nextflow

include { ASSEMBLER_PE } from './../modules/assemblers/main' 
include { CSVTK_UNIQ } from './../modules/csvtk/main'

workflow RUN_ASSEMBLE {   

    take:
        input
        
    main:
    // reads = input.filter{ cfg,files -> cfg.modules.contains("assemble") && cfg.pe_reads == true }
    // contigs = input.filter{ cfg,files -> !cfg.modules.contains("assemble") && cfg.asm == true }
    pe_reads = input.filter{ cfg,files -> cfg.input_type == 'pe_reads' }
    
    ASSEMBLER_PE ( pe_reads )
    contigs = ASSEMBLER_PE.out.contigs
    contigs = contigs.map { cfg, files -> tuple(cfg + [input_type:"asm"], files) }
    // log = ASSEMBLER_PE.out.log
    // flye = reads.filter{ cfg,files -> cfg.assembler == 'flye' }
    versions = ASSEMBLER_PE.out.version.map { cfg, file -> file }.collect()
                                        .map { files -> tuple("version_assembler", files) }
    CSVTK_UNIQ ( versions )
    
    
    // FLYE ( flye )
    emit:
        contigs
        versions = CSVTK_UNIQ.out.collated
        
        
}