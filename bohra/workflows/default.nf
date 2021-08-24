#!/usr/bin/env nextflow

include { SNIPPY } from './../modules/snippy/main' 
include { SNIPPY_CORE } from './../modules/snippy_core/main' 
include { SNP_DISTS } from './../modules/snp_dists/main' 
include { SNIPPY_QC;ADD_HEADER_MLST;MOBSUITE_WRANGLE;COLLATE_MOBSUITE } from './../modules/collation/main' 
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
        SNP_DISTS ( SNIPPY_CORE.out.core_full_aln )     
        
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
        GUBBINS ( core_full )
    emit:
        core_aln = GUBBINS.out.gubbins
}

workflow RUN_IQTREE {
    take:
        core_aln
        reference
    main:
        IQTREE ( core_aln,reference )
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
        ADD_HEADER_MLST ( MLST.out.mlst)
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
        // add the result file emmitted by new rule
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

