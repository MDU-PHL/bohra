#!/usr/bin/env nextflow

include { SNIPPY } from './../modules/snippy/main' 
include { SNIPPY_CORE } from './../modules/snippy_core/main' 
include { SNP_DISTS } from './../modules/snp_dists/main' 
include { SNIPPY_QC } from './../modules/collation/main' 
include { GUBBINS } from './../modules/gubbins/main' 
include { SNIPPY_CLEAN } from './../modules/snippy_clean/main' 
include { IQTREE } from './../modules/iqtree/main' 
include { CSVTK_CONCAT } from './../modules/csvtk/main'


workflow RUN_SNPS {

    take:
        reads
        reference
    main:
        SNIPPY ( reads.combine( reference ) )  
        SNIPPY_QC ( SNIPPY.out.aln )
        alns =  SNIPPY.out.aln.map { cfg, aln -> aln.getParent() }.collect()
        SNIPPY_CORE ( alns, reference )  
        
        core_aln =  SNIPPY_CORE.out.core_aln

        if ( params.gubbins ){
            SNIPPY_CLEAN ( SNIPPY_CORE.out.core_full_aln )
            GUBBINS ( SNIPPY_CLEAN.out.cleaned )
            core_aln = GUBBINS.out.gubbins
        }
        
        SNP_DISTS ( core_aln )

        if (params.run_iqtree ){
            IQTREE ( core_aln, SNIPPY_CORE.out.core_full_aln)
            tree = IQTREE.out.newick
        } else {
            tree = Channel.empty().ifEmpty('EmptyFile')
        }
        stats = SNIPPY_QC.out.snippy_qc.map { cfg, core_stats -> core_stats }.collect().map { files -> tuple("core_genome", files)}
        all_core_stats = CSVTK_CONCAT ( stats )
    emit:
        
        dists = SNP_DISTS.out.distances
        tree = tree
        core_stats = all_core_stats
}

