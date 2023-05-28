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
        SNP_DISTS ( SNIPPY_CORE.out.core_aln )

        core_aln =  SNIPPY_CORE.out.core_aln

        if ( params.gubbins ){
            SNIPPY_CLEAN ( SNIPPY_CORE.out.core_full )
            GUBBINS ( SNIPPY_CLEAN.out.cleaned )
            core_aln = GUBBINS.out.gubbins
        }
        if (params.run_iqtree ){
            IQTREE ( core_aln, SNIPPY_CORE.out.core_full_aln)
            tree = IQTREE.out.newick
        } else {
            tree = Channel.empty().ifEmpty('EmptyFile')
        }
        stats = SNIPPY_QC.out.snippy_qc.map { cfg, core_stats -> core_stats }.collect().map { files -> tuple("core_genome", files)}
        all_core_stats = CSVTK_CONCAT ( stats )
    emit:
        // aln = SNIPPY.out.aln
        // qual = SNIPPY_QC.out.snippy_qc
        // core_aln = SNIPPY_CORE.out.core_aln
        // core_full_aln = SNIPPY_CORE.out.core_full_aln
        // core_vcf = SNIPPY_CORE.out.core_vcf
        // core_stats = SNIPPY_CORE.out.core_stats
        dists = SNP_DISTS.out.distances
        tree = tree
        core_stats = all_core_stats
}

// workflow RUN_SNIPPY {   

//     take:
//         reads
//     main:
//         SNIPPY ( reads )  
//         SNIPPY_QC ( SNIPPY.out.aln )
        
//     emit:
//         aln = SNIPPY.out.aln
//         qual = SNIPPY_QC.out.snippy_qc

// }


// workflow RUN_CORE {   

//     take:
//         alns
        
//     main:
//         SNIPPY_CORE ( alns, reference )  
//         SNP_DISTS ( SNIPPY_CORE.out.core_aln )     
        
//     emit:
//         core_aln = SNIPPY_CORE.out.core_aln
//         core_full_aln = SNIPPY_CORE.out.core_full_aln
//         core_vcf = SNIPPY_CORE.out.core_vcf
//         core_stats = SNIPPY_CORE.out.core_stats
//         dists = SNP_DISTS.out.distances
        
// }

// workflow RUN_IQTREE {
//     take:
//         core_aln
//         core_full_aln
//     main:
//         IQTREE ( core_aln,core_full_aln )
//     emit:
//         newick = IQTREE.out.newick
// }

// workflow RUN_GUBBINS {
//     take:
//         core_full
//     main:
//         SNIPPY_CLEAN ( core_full )
//         GUBBINS ( SNIPPY_CLEAN.out.cleaned )
//     emit:
//         core_aln = GUBBINS.out.gubbins
// }


// workflow CONCAT_CORE_STATS {
//     take:
//         stats
//     main:
//         CSVTK_CONCAT ( stats )
//     emit:
//         collated_core = CSVTK_CONCAT.out.collated

// }

