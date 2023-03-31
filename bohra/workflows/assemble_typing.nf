#!/usr/bin/env nextflow

// include { CSVTK_CONCAT } from './../modules/csvtk/main'
// include { COLLATE_ASM } from './../modules/collation/main'
// include { SHOVILL } from './../modules/shovill/main' 
// include { SPADES } from './../modules/spades/main' 
// include { SKESA } from './../modules/skesa/main' 
// include { SEQKIT_STATS } from './../modules/seqkit/stats/main' addParams( options: [args2: 'contigs'] )
// include { ABRITAMR;COMBINE_AMR } from './../modules/resistome/main' addParams( options: [args2: 4] )
// include { MLST } from './../modules/mlst/main' addParams( options: [args2: 4] )
// include { PROKKA } from './../modules/prokka/main' 

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
        println contigs.view()
        // SEQKIT_STATS ( contigs )
        // ABRITAMR ( contigs )
        // ab = ABRITAMR.out.abritamr_matches.join(ABRITAMR.out.abritamr_partials)
        // COMBINE_AMR( ab )
        // MLST ( contigs )
        // PROKKA ( contigs )
        
    // emit:
    //     contigs
    //     assembly_stats = SEQKIT_STATS.out.stats
    //     resistome = COMBINE_ABRITAMR.out.resistome
    //     virulome = ABRITAMR.out.abritamr_virulence
    //     mlst = MLST.out.mlst
    //     gff = PROKKA.out.gff
    //     prokka_txt = PROKKA.out.prokka_txt
}

// workflow CONCAT_MLST {
//     take:
//         mlsts
//     main:
//         CSVTK_CONCAT ( mlsts )
//     emit:
//         collated_mlst = CSVTK_CONCAT.out.collated

// }


// workflow CONCAT_RESISTOMES {
//     take:
//         resistomes
//     main:
//         CSVTK_CONCAT ( resistomes )
//     emit:
//         collated_resistomes = CSVTK_CONCAT.out.collated

// }


// workflow CONCAT_ASM {
//     take:
//         asm
//     main:
//         CSVTK_CONCAT ( asm )
//     emit:
//         collated_assembly = CSVTK_CONCAT.out.collated

// }

// workflow CONCAT_VIRULOME {
//     take:
//         virulome
//     main:
//         CSVTK_CONCAT ( virulome )
//     emit:
//         collated_virulome = CSVTK_CONCAT.out.collated

// }


// workflow COLLATE_ASM_PROKKA {
//     take:
//         asm
//     main:
//         COLLATE_ASM ( asm )
//     emit:
//         collated_asm = COLLATE_ASM.out.assembly

// }