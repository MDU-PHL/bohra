#!/usr/bin/env nextflow

include { SNIPPY } from './../modules/snippy/main' 
include { SNIPPY_CORE } from './../modules/snippy_core/main' 
include { SNP_DISTS } from './../modules/snp_dists/main' 
include { SNIPPY_QC } from './../modules/collation/main' 
include { GUBBINS } from './../modules/gubbins/main' 
include { SNIPPY_CLEAN } from './../modules/snippy_clean/main' 


workflow RUN_SNIPPY {   

    take:
        reads
    main:
        SNIPPY ( reads )  
        SNIPPY_QC ( SNIPPY.out.aln )
    emit:
        aln = SNIPPY.out.aln
        qual = SNIPPY_QC.out.snippy_qc

}


workflow RUN_CORE {   

    take:
        alns
        reference
    main:
        SNIPPY_CORE ( alns, reference )  
        SNP_DISTS ( SNIPPY_CORE.out.core_aln )     
        
    emit:
        core_aln = SNIPPY_CORE.out.core_aln
        core_full_aln = SNIPPY_CORE.out.core_full_aln
        core_vcf = SNIPPY_CORE.out.core_vcf
        core_stats = SNIPPY_CORE.out.core_stats
        dists = SNP_DISTS.out.distances
        
}

workflow RUN_GUBBINS {
    take:
        core_full
    main:
        SNIPPY_CLEAN ( core_full )
        GUBBINS ( SNIPPY_CLEAN.out.cleaned )
    emit:
        core_aln = GUBBINS.out.gubbins
}


workflow CONCAT_CORE_STATS {
    take:
        stats
    main:
        CSVTK_CONCAT ( stats )
    emit:
        collated_core = CSVTK_CONCAT.out.collated

}

