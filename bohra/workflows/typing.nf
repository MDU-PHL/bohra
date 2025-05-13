#!/usr/bin/env nextflow

// include { ABRITAMR;COMBINE_AMR } from './../modules/resistome/main' addParams( options: [args2: 4] )
// include { MLST } from './../modules/mlst/main' addParams( options: [args2: 4] )
// include { MOBSUITE } from './../modules/mobsuite/main' 
// include { ADD_HEADER_MLST;COLLATE_ABRITMAR;COLLATE_MOBSUITE;MOBSUITE_WRANGLE } from './../modules/collation/main' 
include { SEROTYPES } from './../subworkflows/serotypes'
include { RUN_MLST } from './../subworkflows/mlst'
include {RUN_MOBSUITE} from './../subworkflows/plasmids'
include { RUN_ABRITAMR;CONCAT_INFERRED;CONCAT_PLASMID;CONCAT_REPORTABLE;CONCAT_RESISTOMES;CONCAT_VIRULENCE } from './../subworkflows/abritamr'
// include {CSVTK_CONCAT } from './../modules/csvtk/main'


workflow RUN_TYPING {

    take:
        asm
        reads
    main:
        SEROTYPES ( asm, reads )
        serotypes = SEROTYPES.out.collated_typers
        RUN_MLST ( asm )
        mlst = RUN_MLST.out.collated_mlst
        RUN_MOBSUITE ( asm )
        RUN_ABRITAMR ( asm, RUN_MOBSUITE.out.contig_report )


        resistomes= RUN_ABRITAMR.out.resistome.map{ cfg, file ->file}.collect()
        resistomes = resistomes.map { files -> tuple("resistome", files) }
        virulences = RUN_ABRITAMR.out.abritamr_virulence.map{ cfg, file -> file}.collect()
        virulences = virulences.map { files -> tuple("virulence", files) }
        plasmids = RUN_ABRITAMR.out.plasmid.map{ cfg, file ->file}.collect()
        plasmids = plasmids.map { files -> tuple("plasmid", files) }
        inferreds = RUN_ABRITAMR.out.abritamr_infer.map{ cfg, file -> file}.collect()
        inferreds = inferreds.map { files -> tuple("inferred_antibiogram", files) }
        reportables = RUN_ABRITAMR.out.abritamr_reportable.map{ cfg, file -> file}.collect()
        reportables = reportables.map { files -> tuple("reportable_amr_genes", files) }
        // println plasmids.view()
        resistome = CONCAT_RESISTOMES ( resistomes )
        virulome = CONCAT_VIRULENCE ( virulences )
        plasmid = CONCAT_PLASMID ( plasmids )
        inferred = CONCAT_INFERRED ( inferreds )
        reportable = CONCAT_REPORTABLE ( reportables )
        // species_stats.map { files -> tuple("speciation", files) }
    emit:
        resistome = resistome
        virulome = virulome
        plasmid = plasmid
        inferred = inferred
        reportable = reportable
        serotypes = serotypes
        mlst = mlst

}


// workflow CONCAT_TYPER {

//     take:
//         typers
//     main:
//         CSVTK_CONCAT ( typers )
//     emit:
//         collated_typers = CSVTK_CONCAT.out.collated

// }

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
//         COLLATE_ABRITMAR ( resistomes )
//     emit:
//         collated_resistomes = COLLATE_ABRITMAR.out.collated

// }

// workflow CONCAT_VIRULOMES {
//     take:
//         virulomes
//     main:
//         COLLATE_ABRITMAR ( virulomes )
//     emit:
//         collated_virulomes = COLLATE_ABRITMAR.out.collated

// }


// workflow CONCAT_PLASMID {
//     take:
//         plasmids
//     main:
//         COLLATE_MOBSUITE ( plasmids )
//     emit:
//         collated_plasmids = COLLATE_MOBSUITE.out.collated_plasmid

// }


// workflow BASIC_TYPING {
    
//     take:
//         contigs
//     main:
//         ABRITAMR ( contigs )
//         ab = ABRITAMR.out.abritamr_matches.join(ABRITAMR.out.abritamr_partials)
//         COMBINE_AMR( ab )
//         MLST ( contigs )
//         ADD_HEADER_MLST ( MLST.out.json)
//         MOBSUITE ( contigs )
//         MOBSUITE_WRANGLE ( MOBSUITE.out.mobs )
        
//     emit:
//         resistome = CONCAT_RESISTOMES ( COMBINE_AMR.out.resistome.map { cfg, resistome -> resistome }.collect().map { files -> tuple("resistome", files)})
//         virulome = CONCAT_VIRULOMES ( ABRITAMR.out.abritamr_virulence.map { cfg, virulome -> virulome }.collect().map { files -> tuple("virulome", files)})
//         mlst = CONCAT_MLST ( ADD_HEADER_MLST.out.mlst.map { cfg, mlst -> mlst }.collect().map { files -> tuple("mlst", files)})
//         plasmid = CONCAT_PLASMID ( MOBSUITE_WRANGLE.out.plasmid.map { cfg, plasmid -> plasmid }.collect().map { files -> tuple("plasmid", files)})
        

// }


