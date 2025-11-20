#!/usr/bin/env nextflow

include { ASSEMBLER_PE } from './../modules/assemblers/main' 
include { CSVTK_UNIQ } from './../modules/csvtk/main'
include { INSERTIQR } from './../modules/insertiqr/main'
include { CONCAT_FILES } from './../modules/utils/main'


workflow RUN_ASSEMBLE {   

    take:
        input
        
    main:
    pe_reads = input.filter{ cfg,files -> cfg.input_type == 'pe_reads' }
    // println pe_reads.view()
    ASSEMBLER_PE ( pe_reads )
    contigs = ASSEMBLER_PE.out.contigs
    contigs = contigs.map { cfg, files -> tuple(cfg + [input_type:"asm"], files) }
    ctg = contigs.map { cfg, file -> tuple([id:cfg.id], file) }
    per = pe_reads.map { cfg, files -> tuple([id:cfg.id], files)}
    forinsert = ctg.join( per )
    // println forinsert.view()
                                        // .map { id, cfg_asm, files, cfg_pe, reads -> tuple(cfg_asm + [pe_reads:reads], files) }
    INSERTIQR ( forinsert )
    // log = ASSEMBLER_PE.out.log
    // flye = reads.filter{ cfg,files -> cfg.assembler == 'flye' }
    versions = ASSEMBLER_PE.out.version.map { cfg, file -> file }.collect()
                                        .map { files -> tuple("version_assembler", files) }
    CSVTK_UNIQ ( versions )
    
    insertiqr = INSERTIQR.out.stats.map{ cfg, file -> file}.collect()
    insertiqr = insertiqr   .map { files -> tuple("insertiqr", files) }
    CONCAT_FILES ( insertiqr )
    // FLYE ( flye )
    emit:
        contigs
        insertiqr = CONCAT_FILES.out.collated
        versions = CSVTK_UNIQ.out.collated
        
        
}