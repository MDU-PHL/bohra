#!/usr/bin/env nextflow
include { LISSERO } from './../modules/lissero/main'
include { MENINGOTYPE } from './../modules/meningotype/main'
include { STYPE } from './../modules/stype/main'
include { NGMASTER } from './../modules/ngmaster/main'
include { KLEBORATE } from './../modules/kleborate/main'
include { ECTYPER } from './../modules/ectyper/main'
include { EMMTYPER } from './../modules/emmtyper/main'
include {CSVTK_CONCAT;CSVTK_UNIQ } from './../modules/csvtk/main'

workflow SEROTYPES {
    
    take:
        asm
        reads
        
    main:
        
        listeria = asm.filter { cfg, contigs -> cfg.species == 'Listeria monocytogenes'}
        nmen = asm.filter { cfg, contigs -> cfg.species == 'Neisseria meningitidis'}
        ngono = asm.filter { cfg, contigs -> cfg.species == 'Neisseria gonorrhoeae'}
        salmonella = asm.filter { cfg, contigs -> cfg.species =~ 'Salmonella'}
        klebs = asm.filter { cfg, contigs -> cfg.species =~ 'Klebsiella'}
        ecoli = asm.filter { cfg, contigs -> cfg.species == 'Escherichia coli'}
        igas = asm.filter { cfg, contigs -> cfg.species == 'Streptococcus pyogenes'}
        // add in shigella
        LISSERO ( listeria )
        lissero_typers = LISSERO.out.typer.map {cfg, typer -> typer }.collect()
        lissero_version = LISSERO.out.version.map {cfg, version -> version }.collect()
        STYPE ( salmonella )
        salmo_typers = STYPE.out.typer.map {cfg, typer -> typer }.collect()
        salmo_version = STYPE.out.version.map {cfg, version -> version }.collect()
        MENINGOTYPE ( nmen )
        nmen_typers = MENINGOTYPE.out.typer.map {cfg, typer -> typer }.collect()
        nmen_version = MENINGOTYPE.out.version.map {cfg, version -> version }.collect()
        NGMASTER ( ngono )
        ngono_typers = NGMASTER.out.typer.map {cfg, typer -> typer }.collect()
        ngono_version = NGMASTER.out.version.map {cfg, version -> version }.collect()
        KLEBORATE ( klebs )
        klebs_typers = KLEBORATE.out.typer.map {cfg, typer -> typer }.collect()
        klebs_version = KLEBORATE.out.version.map {cfg, version -> version }.collect()
        ECTYPER ( ecoli )
        ecoli_typers = ECTYPER.out.typer.map {cfg, typer -> typer }.collect()
        ecoli_version = ECTYPER.out.version.map {cfg, version -> version }.collect()
        EMMTYPER ( igas )
        emm_typers = EMMTYPER.out.typer.map {cfg, typer -> typer }.collect()
        emm_version = EMMTYPER.out.version.map {cfg, version -> version }.collect()
        typers = lissero_typers.concat ( salmo_typers, nmen_typers, ngono_typers, klebs_typers, ecoli_typers,emm_typers ).map { files -> tuple("typer", files)}
        versions = lissero_version.concat ( salmo_version, nmen_version, ngono_version, klebs_version, ecoli_version, emm_version ).map { files -> tuple("version_serotypes", files)}
        CSVTK_CONCAT ( typers )
        CSVTK_UNIQ ( versions )
        collated_versions = CSVTK_UNIQ.out.collated
        collated_typers = CSVTK_CONCAT.out.collated
    emit:
        collated_typers
        collated_versions

}
