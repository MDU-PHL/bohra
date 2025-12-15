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
include { FILTER_CORE  } from './../modules/check_core/main'
include { CHECK_CORE } from './../modules/check_core/main'
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
        snippy_qc = SNIPPY_QC.out.snippy_qc.map { cfg, snippy_qc -> snippy_qc }.collect().map { it.join(' ') }
        // println snippy_qc.view()
        // println snippy_qc.channelType()
        CHECK_CORE ( snippy_qc )
        all_core_stats = CHECK_CORE.out.stats
    //     // need to come up with a way to add a flag to a aln to prevent inclusion.... 
        FILTER_CORE ( SNIPPY.out.aln.combine( all_core_stats ))

        
        alns_qc = FILTER_CORE.out.aln_filter.join ( SNIPPY.out.aln )
        alns_qc = alns_qc.map { cfg, val, aln -> tuple( cfg + [filter:val.trim()]), aln}
        
        
        if ( params.ignore_warnings == true) {
            println "Including all sequences in core genome alignment regardless of QC metrics."
            alns = alns_qc.map { cfg, aln -> aln.getParent() }.collect()
        } else {
            println "Excluding outlier sequences from core genome alignment based on QC metrics."
            
            alns = alns_qc.filter { cfg, aln -> cfg.filter == "include"}.map { cfg, aln -> aln.getParent() }.collect()
            
        }
               

        SNIPPY_CORE ( alns, reference )  
        CORE_SNP_FILTER ( SNIPPY_CORE.out.core_full_aln )
        core_aln =  CORE_SNP_FILTER.out.aln
        core_full_aln = SNIPPY_CORE.out.core_full_aln
        core_vcf = SNIPPY_CORE.out.core_vcf
        SNIPPY_CLEAN ( core_full_aln )
        cleaned_aln = SNIPPY_CLEAN.out.cleaned
        
        if ( params.gubbins ){
            GUBBINS ( cleaned_aln )
            core_aln = GUBBINS.out.gubbins
        }
        
        SNP_DISTS ( core_aln )
        if ( params.cluster ) {
            SNP_CLUSTER ( SNP_DISTS.out.matrix )
        } else {
            SNP_CLUSTER.out.clusters = Chanel.empty().ifEmpty("not_available")
        }
        
    emit:
        
        dists = SNP_DISTS.out.matrix
        aln = core_aln
        core_full_aln = core_full_aln
        core_vcf = core_vcf
        // cleaned_aln = cleaned_aln
        stats = all_core_stats
        clusters = SNP_CLUSTER.out.clusters
        version = CSVTK_UNIQ.out.collated

}