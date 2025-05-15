#!/usr/bin/env nextflow

include { SEQTK } from './../modules/seqtk/main' 
include { COLLATE_STATS_ISOLATE;COLLATE_ASM } from './../modules/collation/main'
include { SEQKIT_STATS } from './../modules/seqkit/stats/main' 
include { SEQKIT_GC } from './../modules/seqkit/fx2tab/main' 
include { KMC } from './../modules/kmc/main' 
include { PROKKA } from './../modules/prokka/main'
include { CSVTK_CONCAT;CSVTK_UNIQ } from './../modules/csvtk/main'

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
        versions_seqkit = SEQKIT_STATS.out.version.map { cfg, file -> file }.collect()
                                         .map { files -> tuple("version_seqkit", files) }
        VERSION_SEQKIT_READS ( versions_seqkit )
        versions_kmc = KMC.out.version.map { cfg, file -> file }.collect()
                                         .map { files -> tuple("version_kmc", files) }
        VERSION_KMC ( versions_kmc )
        CSVTK_CONCAT ( seq_stats )
    emit:
        read_stats  = CSVTK_CONCAT.out.collated
        version_seqkit_reads = VERSION_SEQKIT_READS.out.version
        version_kmc = VERSION_KMC.out.version
        

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
        versions_prokka = PROKKA.out.version.map { cfg, file -> file }.collect()
                                         .map { files -> tuple("version_prokka", files) }
        versions_seqkit = SEQKIT_STATS.out.version.map { cfg, file -> file }.collect()
                                         .map { files -> tuple("version_seqkit", files) }
        VERSION_PROKKA ( versions_prokka )
        VERSION_SEQKIT_ASM ( versions_seqkit )

        CSVTK_CONCAT ( asm_stats )
    emit:
        assembly_stats = CSVTK_CONCAT.out.collated
        version_prokka = VERSION_PROKKA.out.version
        version_seqkit_asm = VERSION_SEQKIT_ASM.out.version
        
        
}

workflow VERSION_PROKKA {
    take:
        prokka
    main:
        
        CSVTK_UNIQ ( prokka )
    emit:
        version = CSVTK_UNIQ.out.collated
}

workflow VERSION_KMC {
    take:
        kmc
    main:
        
        CSVTK_UNIQ ( kmc )
    emit:
        version = CSVTK_UNIQ.out.collated
}

workflow VERSION_SEQKIT_READS {
    take:
        seqkit
    main:
        
        CSVTK_UNIQ ( seqkit )
    emit:
        version = CSVTK_UNIQ.out.collated
}


workflow VERSION_SEQKIT_ASM {
    take:
        seqkit
    main:
        
        CSVTK_UNIQ ( seqkit )
    emit:
        version = CSVTK_UNIQ.out.collated
}