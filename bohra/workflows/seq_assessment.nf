#!/usr/bin/env nextflow

include { SEQTK } from './../modules/seqtk/main' 
include { COLLATE_STATS_ISOLATE;COLLATE_ASM } from './../modules/collation/main'
include { SEQKIT_STATS } from './../modules/seqkit/stats/main' 
include { SEQKIT_GC } from './../modules/seqkit/fx2tab/main' 
include { KMC } from './../modules/kmc/main' 
include { PROKKA } from './../modules/prokka/main'
include { CSVTK_CONCAT } from './../modules/csvtk/main'

workflow READ_ANALYSIS {   

    take:
        reads_pe
        
        
    main:
        
        SEQKIT_STATS ( reads_pe )
        SEQKIT_GC ( reads_pe )
        KMC ( reads_pe )
        COMBD = SEQKIT_STATS.out.stats.join( SEQKIT_GC.out.stats )
        COMBD = COMBD.join( KMC.out.genome_size )
        COLLATE_STATS_ISOLATE ( COMBD )
        seq_stats = COLLATE_STATS_ISOLATE.out.read_assessment.map { cfg, seq -> seq }.collect()
        seq_stats = seq_stats.map { files -> tuple("read_assessment", files) }
        CSVTK_CONCAT ( seq_stats )
    emit:
        read_stats  = CSVTK_CONCAT.out.collated
        

}


workflow ASSEMBLY_ANALYSIS {   

    take:
        contigs
        
    main:
        SEQKIT_GC ( contigs )
        SEQKIT_STATS ( contigs )
        PROKKA ( contigs )
        APS = PROKKA.out.prokka_txt.join( SEQKIT_STATS.out.stats )
        COLLATE_ASM ( APS )
        asm_stats = COLLATE_ASM.out.assembly.map { cfg, asm -> asm }.collect()
        asm_stats = asm_stats.map { files -> tuple("assembly_assesment", files) }
        
        CSVTK_CONCAT ( asm_stats )
    emit:
        assembly_stats = CSVTK_CONCAT.out.collated
        
        
}