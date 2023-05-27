#!/usr/bin/env nextflow

include { ABRITAMR;COMBINE_AMR } from './../modules/resistome/main' addParams( options: [args2: 4] )
include { MLST } from './../modules/mlst/main' addParams( options: [args2: 4] )
// include { PROKKA } from './../modules/prokka/main' 
include { ADD_HEADER_MLST } from './../modules/collation/main' 
include { LISSERO } from './../modules/lissero/main'
include { MENINGOTYPE } from './../modules/meningotype/main'
include { STYPE } from './../modules/stype/main'
include {CSVTK_CONCAT } from './../modules/csvtk/main'
workflow BASIC {
    
    take:
        contigs
    main:
        ABRITAMR ( contigs )
        ab = ABRITAMR.out.abritamr_matches.join(ABRITAMR.out.abritamr_partials)
        COMBINE_AMR( ab )
        MLST ( contigs )
        ADD_HEADER_MLST ( MLST.out.json)
        
        
    emit:
        resistome = COMBINE_AMR.out.resistome
        virulome = ABRITAMR.out.abritamr_virulence
        mlst = ADD_HEADER_MLST.out.mlst
        

}



workflow SEROTYPES {
    
    take:
        typing_input
        // species
    main:
        
        listeria = typing_input.filter { cfg, contigs, species_obs -> species_obs == 'Listeria monocytogenes'}
        nmen = typing_input.filter { cfg, contigs, species_obs -> species_obs == 'Neisseria meningitidis'}
        ngono = typing_input.filter { cfg, contigs, species_obs -> species_obs == 'Neisseria gonorrhoeae'}
        salmonella = typing_input.filter { cfg, contigs, species_obs -> species_obs =~ 'Salmonella'}
        // listeria = listeria.map { cfg, contigs, species_obs -> cfg,contigs }
        // println listeria.view()
        LISSERO ( listeria )
        lissero_typers = LISSERO.out.typer.map {cfg, typer -> typer }.collect()
        STYPE ( salmonella )
        salmo_typers = STYPE.out.typer.map {cfg, typer -> typer }.collect()
        MENINGOTYPE ( nmen )

        typers = lissero_typers.concat ( salmo_typers, nmen ).map { files -> tuple("typer", files)}
        println typers.view()

    emit:
        typers = typers

}

workflow CONCAT_TYPER {

    take:
        typers
    main:
        CSVTK_CONCAT ( typers )
    emit:
        collated_typers = CSVTK_CONCAT.out.collated

}