#!/usr/bin/env nextflow

include { SNIPPY } from './../modules/snippy/main' 
include { SNIPPY_CORE } from './../modules/snippy_core/main' 
include { SNP_DISTS } from './../modules/snp_dists/main' 
include { SNIPPY_QC } from './../modules/collation/main' 
include { GUBBINS } from './../modules/gubbins/main' 
include { SNIPPY_CLEAN } from './../modules/snippy_clean/main' 
include { CSVTK_CONCAT;CSVTK_UNIQ } from './../modules/csvtk/main'
include { CORE_SNP_FILTER } from './../modules/core_snp_filter/main'
include { SNP_CLUSTER } from './../modules/cluster/main'

workflow RUN_SNPS {

    take:
        reads
        reference
    main:
        
        SNIPPY ( reads.combine( reference ) )  
        versions = SNIPPY.out.version.map { cfg, file -> file }.collect()
                                        .map { files -> tuple("version_snippy", files) }
        CSVTK_UNIQ ( versions )
        SNIPPY_QC ( SNIPPY.out.aln )
        alns =  SNIPPY.out.aln.map { cfg, aln -> aln.getParent() }.collect()
        SNIPPY_CORE ( alns, reference )  
        CORE_SNP_FILTER ( SNIPPY_CORE.out.core_full_aln )
        core_aln =  CORE_SNP_FILTER.out.aln
        core_full_aln = SNIPPY_CORE.out.core_full_aln
        SNIPPY_CLEAN ( core_full_aln )
        cleaned_aln = SNIPPY_CLEAN.out.cleaned
        
        if ( params.gubbins ){
            GUBBINS ( cleaned_aln )
            core_aln = GUBBINS.out.gubbins
        }
        
        SNP_DISTS ( core_aln )
        SNP_CLUSTER ( SNP_DISTS.out.distances )
        stats = SNIPPY_QC.out.snippy_qc.map { cfg, core_stats -> core_stats }.collect().map { files -> tuple("core_genome", files)}
        all_core_stats = CSVTK_CONCAT ( stats )
    emit:
        
        dists = SNP_DISTS.out.distances
        aln = core_aln
        core_full_aln = core_full_aln
        // cleaned_aln = cleaned_aln
        stats = all_core_stats
        clusters = SNP_CLUSTER.out.clusters
        version = CSVTK_UNIQ.out.collated

}