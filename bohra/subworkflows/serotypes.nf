#!/usr/bin/env nextflow
include { LISSERO } from './../modules/lissero/main'
include { MENINGOTYPE } from './../modules/meningotype/main'
include { STYPE } from './../modules/stype/main'
include { NGMASTER } from './../modules/ngmaster/main'
include { KLEBORATE } from './../modules/kleborate/main'
include { ECTYPER } from './../modules/ectyper/main'
include { EMMTYPER } from './../modules/emmtyper/main'
include {CSVTK_CONCAT } from './../modules/csvtk/main'

workflow SEROTYPES {
    
    take:
        typing_input
        
    main:
        
        listeria = typing_input.filter { cfg, contigs -> cfg.species == 'Listeria monocytogenes'}
        nmen = typing_input.filter { cfg, contigs -> cfg.species == 'Neisseria meningitidis'}
        ngono = typing_input.filter { cfg, contigs -> cfg.species == 'Neisseria gonorrhoeae'}
        salmonella = typing_input.filter { cfg, contigs -> cfg.species =~ 'Salmonella'}
        klebs = typing_input.filter { cfg, contigs -> cfg.species =~ 'Klebsiella'}
        ecoli = typing_input.filter { cfg, contigs -> cfg.species == 'Escherichia coli'}
        igas = typing_input.filter { cfg, contigs -> cfg.species == 'Streptococcus pyogenes'}
        // add in shigella
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
        CSVTK_CONCAT ( typers )
        collated_typers = CSVTK_CONCAT.out.collated
    emit:
        collated_typers

}
