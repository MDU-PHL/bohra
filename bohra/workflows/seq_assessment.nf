#!/usr/bin/env nextflow

include { SEQTK } from './../modules/seqtk/main' 
include { COLLATE_STATS_ISOLATE;COLLATE_ASM } from './../modules/collation/main'
include { SEQKIT_STATS } from './../modules/seqkit/stats/main' 
include { SEQKIT_GC } from './../modules/seqkit/fx2tab/main' 
include { KMC } from './../modules/kmc/main' 
include { PROKKA } from './../modules/prokka/main'
include { CSVTK_CONCAT;CSVTK_UNIQ } from './../modules/csvtk/main'
include { BOHRA_VERSION } from './../modules/utils/main'
include { CHECK_FASTQ } from './../modules/check_fastq/main'
include { FASTP } from './../modules/fastp/main'

workflow READ_ANALYSIS {   

    take:
        reads_pe
        
        
    main:
        BOHRA_VERSION (  )
        if (params.modules.contains('trim')){
            FASTP( reads_pe )
            reads_pe = FASTP.out.read

        }
        CHECK_FASTQ ( reads_pe )
        // println reads_pe.view()
        reads_pe_checked = CHECK_FASTQ.out.pe_check.join( reads_pe )
        reads_pe_checked = reads_pe_checked.map { cfg, check, files -> tuple(cfg + [check:check], files) }
        
        // println reads_pe_checked.view()
        SEQKIT_STATS ( reads_pe_checked )
        SEQKIT_GC ( reads_pe_checked )
        KMC ( reads_pe_checked )
        SEQTK ( reads_pe_checked )
        COMBD = SEQKIT_STATS.out.stats.join( SEQKIT_GC.out.stats )
        COMBD = COMBD.join( KMC.out.genome_size )
        COMBD = COMBD.join( SEQTK.out.seqtk_stats )
        COLLATE_STATS_ISOLATE ( COMBD )
        seq_stats = COLLATE_STATS_ISOLATE.out.read_assessment.map { cfg, seq -> seq }.collect() // TODO add in Q score
        seq_stats = seq_stats.map { files -> tuple("read_assessment", files) }
        versions_seqkit = SEQKIT_STATS.out.version.map { cfg, file -> file }.collect()
                                         .map { files -> tuple("version_seqkit", files) }
        VERSION_SEQKIT_READS ( versions_seqkit )
        versions_kmc = KMC.out.version.map { cfg, file -> file }.collect()
                                         .map { files -> tuple("version_kmc", files) }
        VERSION_KMC ( versions_kmc )
        CSVTK_CONCAT ( seq_stats )
        reads_pe_checked = reads_pe_checked.filter { cfg, files -> cfg.check != 'FAIL_READ_FILE_TOO_SMALL' }
                                            .map { cfg,files -> tuple(cfg.findAll {it.key != 'check'}, files ) }
    emit:
        read_stats  = CSVTK_CONCAT.out.collated
        reads_pe = reads_pe_checked
        version_seqkit_reads = VERSION_SEQKIT_READS.out.version
        version_kmc = VERSION_KMC.out.version
        version_bohra = BOHRA_VERSION.out.collated
        

}


workflow ASSEMBLY_ANALYSIS {   

    take:
        contigs
        
    main:
        // println contigs.view()
        BOHRA_VERSION (  )
        SEQKIT_GC ( contigs )
        SEQKIT_STATS ( contigs )
        // println SEQKIT_STATS.out.stats.view()
        PROKKA ( contigs )
        gff = PROKKA.out.gff.filter { cfg, files -> cfg.control != 'control' }
        gff = gff.map { cfg, files -> files }.collect()
        APS = PROKKA.out.prokka_txt.join( SEQKIT_STATS.out.stats )
        COLLATE_ASM ( APS )
        asm_stats = COLLATE_ASM.out.assembly.map { cfg, asm -> asm }.collect()
        asm_stats = asm_stats.map { files -> tuple("assembly_assessment", files) }
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
        version_bohra = BOHRA_VERSION.out.collated
        gff
        
        
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