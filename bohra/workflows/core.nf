#!/usr/bin/env nextflow

include { SNIPPY_CORE } from './../modules/snippy_core/main' 
include { IQTREE } from './../modules/iqtree/main' addParams( options: [args2: params.iqtree_cpus] )

workflow RUN_CORE {   

    take:
        alns
    main:
        SNIPPY_CORE ( alns )       
        IQTREE ( SNIPPY_CORE.out.core_aln )
    emit:
        core_aln = SNIPPY_CORE.out.core_aln
        core_full_aln = SNIPPY_CORE.out.core_full_aln
        core_vcf = SNIPPY_CORE.out.core_vcf
        core_stats = SNIPPY_CORE.out.core_stats
        newick = IQTREE.out.newick


}

