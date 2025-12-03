#!/usr/bin/env nextflow

// include { ABRITAMR;COMBINE_AMR } from './../modules/resistome/main' addParams( options: [args2: 4] )
// include { MLST } from './../modules/mlst/main' addParams( options: [args2: 4] )
// include { MOBSUITE } from './../modules/mobsuite/main' 
// include { ADD_HEADER_MLST;COLLATE_ABRITMAR;COLLATE_MOBSUITE;MOBSUITE_WRANGLE } from './../modules/collation/main' 
include { SEROTYPES } from './../subworkflows/serotypes'
include { RUN_MLST } from './../subworkflows/mlst'
include {RUN_MOBSUITE} from './../subworkflows/plasmids'
include { RUN_ABRITAMR;CONCAT_INFERRED;CONCAT_PLASMID;CONCAT_REPORTABLE;CONCAT_RESISTOMES;CONCAT_VIRULENCE } from './../subworkflows/abritamr'
include { RUN_TBTAMR } from './../workflows/tbtamr'
// include {CSVTK_CONCAT } from './../modules/csvtk/main'


workflow RUN_TYPING {

    take:
        asm
        reads
    main:
        
        RUN_MLST ( asm )
        mlst = RUN_MLST.out.collated_mlst
        SEROTYPES ( asm, reads )
        serotypes = SEROTYPES.out.collated_typers
        
        RUN_MOBSUITE ( asm )
        RUN_ABRITAMR ( asm, RUN_MOBSUITE.out.contig_report )
        versions = RUN_ABRITAMR.out.version.concat( RUN_MOBSUITE.out.version, RUN_MLST.out.version, SEROTYPES.out.collated_versions )
        resistomes= RUN_ABRITAMR.out.resistome.map{ cfg, file ->file}.collect()
        resistomes = resistomes.map { files -> tuple("resistome", files) }
        virulences = RUN_ABRITAMR.out.abritamr_virulence.map{ cfg, file -> file}.collect()
        virulences = virulences.map { files -> tuple("virulence", files) }
        plasmids = RUN_ABRITAMR.out.plasmid.map{ cfg, file ->file}.collect()
        plasmids = plasmids.map { files -> tuple("plasmid", files) }
        inferreds = RUN_ABRITAMR.out.abritamr_infer.map{ cfg, file -> file}.collect()
        inferreds = inferreds.map { files -> tuple("inferred_antibiogram", files) }
        reportables = RUN_ABRITAMR.out.abritamr_reportable.map{ cfg, file -> file}.collect()
        reportables = reportables.map { files -> tuple("reportable_amr_mechanisms", files) }
        // println plasmids.view()
        resistome = CONCAT_RESISTOMES ( resistomes )
        virulome = CONCAT_VIRULENCE ( virulences )
        plasmid = CONCAT_PLASMID ( plasmids )
        inferred = CONCAT_INFERRED ( inferreds )
        reportable = CONCAT_REPORTABLE ( reportables )
        
    emit:
        resistome = resistome
        virulome = virulome
        plasmid = plasmid
        inferred = inferred
        reportable = reportable
        serotypes = serotypes
        mlst = mlst
        versions = versions
        

}
