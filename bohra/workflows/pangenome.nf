#!/usr/bin/env nextflow

include { PANAROO } from './../modules/panaroo/main' 
include { CLASSIFY_PANGENOME } from './../modules/classify_pangenome/main'
include { EXTRACT_GROUPS } from './../modules/extract_groups/main'
workflow RUN_PANAROO {   

    take:
        gff
        group
    main:
        PANAROO ( gff )     
        pangenome_rtab = PANAROO.out.pangenome_rtab
        classification = Channel.empty().ifEmpty( 'no_results' )
        groups = Channel.empty().ifEmpty( 'no_results' )
        if (params.modules.contains("snippy") || (params.modules.contains("ska")) ){
            EXTRACT_GROUPS ( group )
            groups = EXTRACT_GROUPS.out.pangenome_groups
            CLASSIFY_PANGENOME ( pangenome_rtab, groups )
            classification = CLASSIFY_PANGENOME.out.pangenome_classification
        } 
        
    emit:
        classification = classification
        pangenome_rtab = pangenome_rtab
        groups = groups
        roary = PANAROO.out.pangenome_summary
        version = PANAROO.out.version

}