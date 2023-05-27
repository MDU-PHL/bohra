#!/usr/bin/env nextflow

include { ABRITAMR;COMBINE_AMR } from './../modules/resistome/main' addParams( options: [args2: 4] )
include { MLST } from './../modules/mlst/main' addParams( options: [args2: 4] )
// include { PROKKA } from './../modules/prokka/main' 
include { ADD_HEADER_MLST } from './../modules/collation/main' 
include { LISSERO } from './../modules/lissero/main'
include { MENINGOTYPE } from './../modules/meningotype/main'
include { STYPE } from './../modules/stype/main'
include { NGMASTER } from './../modules/ngmaster/main'
include { KLEBORATE } from './../modules/kleborate/main'
include { ECTYPER } from './../modules/ectyper/main'
include { EMMTYPER } from './../modules/emmtyper/main'
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
        
    main:
        
        listeria = typing_input.filter { cfg, contigs, species_obs -> species_obs == 'Listeria monocytogenes'}
        nmen = typing_input.filter { cfg, contigs, species_obs -> species_obs == 'Neisseria meningitidis'}
        ngono = typing_input.filter { cfg, contigs, species_obs -> species_obs == 'Neisseria gonorrhoeae'}
        salmonella = typing_input.filter { cfg, contigs, species_obs -> species_obs =~ 'Salmonella'}
        klebs = typing_input.filter { cfg, contigs, species_obs -> species_obs =~ 'Klebsiella'}
        ecoli = typing_input.filter { cfg, contigs, species_obs -> species_obs == 'Escherichia coli'}
        igas = typing_input.filter { cfg, contigs, species_obs -> species_obs == 'Streptococcus pyogenes'}
        
        LISSERO ( listeria )
        lissero_typers = LISSERO.out.typer.map {cfg, typer -> typer }.collect()
        STYPE ( salmonella )
        salmo_typers = STYPE.out.typer.map {cfg, typer -> typer }.collect()
        MENINGOTYPE ( nmen )
        nmen_typers = MENINGOTYPE.out.typer.map {cfg, typer -> typer }.collect()
        NGMASTER ( ngono )
        ngono_typers = NGMASTER.out.typer.map {cfg, typer -> typer }.collect()
        KLEBORATE ( klebs )
        klebs_typers = KLEBORATE.out.typer.map {cfg, typer -> typer }.collect()
        ECTYPER ( ecoli )
        ecoli_typers = ECTYPER.out.typer.map {cfg, typer -> typer }.collect()
        EMMTYPER ( igas )
        emm_typers = EMMTYPER.out.typer.map {cfg, typer -> typer }.collect()
        typers = lissero_typers.concat ( salmo_typers, nmen_typers, ngono_typers, klebs_typers, ecoli_typers,emm_typers ).map { files -> tuple("typer", files)}

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

