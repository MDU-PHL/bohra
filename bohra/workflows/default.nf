#!/usr/bin/env nextflow

include { SNIPPY } from './../modules/snippy/main' 
include { SNIPPY_CORE } from './../modules/snippy_core/main' 
include { SNIPPY_QC } from './../modules/collation/main' 
include { CSVTK_CONCAT } from './../modules/csvtk/main'
include { COLLATE_ASM } from './../modules/collation/main'
include { SHOVILL } from './../modules/shovill/main' 
include { SPADES } from './../modules/spades/main' 
include { SKESA } from './../modules/skesa/main' 
include { SEQKIT_STATS } from './../modules/seqkit/stats/main' addParams( options: [args2: 'contigs'] )
include { ABRITAMR;COMBINE_AMR } from './../modules/resistome/main' addParams( options: [args2: 4] )
include { MLST } from './../modules/mlst/main' addParams( options: [args2: 4] )
include { PROKKA } from './../modules/prokka/main' 
include { IQTREE } from './../modules/iqtree/main' 


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
    main:
        SNIPPY_CORE ( alns )       
        
    emit:
        core_aln = SNIPPY_CORE.out.core_aln
        core_full_aln = SNIPPY_CORE.out.core_full_aln
        core_vcf = SNIPPY_CORE.out.core_vcf
        core_stats = SNIPPY_CORE.out.core_stats
        
}

workflow RUN_IQTREE {
    take:
        core_aln
    main:
        IQTREE ( core_aln )
    emit:
        newick = IQTREE.out.newick
}
workflow RUN_ASSEMBLE {   

    take:
        reads
    main:
    if ( params.assembler == 'shovill'){
        SHOVILL ( reads )   
        contigs = SHOVILL.out.contigs    
        } 
        else if ( params.assembler == 'spades' ){
        SPADES ( reads )
        contigs = SPADES.out.contigs
        } else if (params.assembler == 'skesa' ) {
        SKESA ( reads )
        contigs = SKESA.out.contigs
        }
        SEQKIT_STATS ( contigs )
        ABRITAMR ( contigs )
        ab = ABRITAMR.out.abritamr_matches.join( ABRITAMR.out.abritamr_partials )
        COMBINE_AMR( ab )
        MLST ( contigs )
        PROKKA ( contigs )
        // println SEQKIT_STATS.out.stats.view()
    emit:
        contigs
        assembly_stats = SEQKIT_STATS.out.stats
        resistome = COMBINE_AMR.out.resistome
        virulome = ABRITAMR.out.abritamr_virulence
        mlst = MLST.out.mlst
        gff = PROKKA.out.gff
        prokka_txt = PROKKA.out.prokka_txt
}

workflow CONCAT_MLST {
    take:
        mlsts
    main:
        CSVTK_CONCAT ( mlsts )
    emit:
        collated_mlst = CSVTK_CONCAT.out.collated

}


workflow CONCAT_RESISTOMES {
    take:
        resistomes
    main:
        CSVTK_CONCAT ( resistomes )
    emit:
        collated_resistomes = CSVTK_CONCAT.out.collated

}

workflow CONCAT_VIRULOMES {
    take:
        virulomes
    main:
        CSVTK_CONCAT ( virulomes )
    emit:
        collated_virulomes = CSVTK_CONCAT.out.collated

}


workflow CONCAT_ASM {
    take:
        asm
    main:
        CSVTK_CONCAT ( asm )
    emit:
        collated_assembly = CSVTK_CONCAT.out.collated

}


workflow COLLATE_ASM_PROKKA {
    take:
        asm
    main:
        COLLATE_ASM ( asm )
    emit:
        collated_asm = COLLATE_ASM.out.assembly

}
