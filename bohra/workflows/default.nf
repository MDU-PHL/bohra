#!/usr/bin/env nextflow

include { SNIPPY } from './../modules/snippy/main' 
include { SNIPPY_CORE } from './../modules/snippy_core/main' 
include { SNP_DISTS } from './../modules/snp_dists/main' 
include { SNIPPY_QC;ADD_HEADER_MLST;MOBSUITE_WRANGLE;COLLATE_MOBSUITE;COLLATE_ABRITMAR } from './../modules/collation/main' 
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
include { MOBSUITE } from './../modules/mobsuite/main' 
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

workflow RUN_IQTREE {
    take:
        core_aln
        core_full_aln
    main:
        IQTREE ( core_aln,core_full_aln )
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
        ADD_HEADER_MLST ( MLST.out.json)
        PROKKA ( contigs )
        // add new processes that generate outputs
        MOBSUITE ( contigs )
        MOBSUITE_WRANGLE ( MOBSUITE.out.mobs )
        // println SEQKIT_STATS.out.stats.view()
    emit:
        contigs
        assembly_stats = SEQKIT_STATS.out.stats
        resistome = COMBINE_AMR.out.resistome
        virulome = ABRITAMR.out.abritamr_virulence
        mlst = ADD_HEADER_MLST.out.mlst
        gff = PROKKA.out.gff
        prokka_txt = PROKKA.out.prokka_txt
        plasmid = MOBSUITE_WRANGLE.out.plasmid
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
        COLLATE_ABRITMAR ( resistomes )
    emit:
        collated_resistomes = COLLATE_ABRITMAR.out.collated

}

workflow CONCAT_VIRULOMES {
    take:
        virulomes
    main:
        COLLATE_ABRITMAR ( virulomes )
    emit:
        collated_virulomes = COLLATE_ABRITMAR.out.collated

}

workflow CONCAT_PLASMID {
    take:
        plasmids
    main:
        COLLATE_MOBSUITE ( plasmids )
    emit:
        collated_plasmids = COLLATE_MOBSUITE.out.collated_plasmid

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


workflow CONCAT_STATS {
    take:
        stats
    main:
        CSVTK_CONCAT ( stats )
    emit:
        collated_stats = CSVTK_CONCAT.out.collated

}

workflow CONCAT_CORE_STATS {
    take:
        stats
    main:
        CSVTK_CONCAT ( stats )
    emit:
        collated_core = CSVTK_CONCAT.out.collated

}

