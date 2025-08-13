#!/usr/bin/env nextflow
include { LISSERO } from './../modules/lissero/main'
include { MENINGOTYPE } from './../modules/meningotype/main'
include { STYPE } from './../modules/stype/main'
include { NGMASTER } from './../modules/ngmaster/main'
include { KLEBORATE } from './../modules/kleborate/main'
include { ECTYPER } from './../modules/ectyper/main'
include { EMMTYPER } from './../modules/emmtyper/main'
include { SHIGAPASS } from './../modules/shigapass/main'
include { SONNEITYPE } from './../modules/sonneitype/main'
include {CSVTK_CONCAT;CSVTK_UNIQ } from './../modules/csvtk/main'
include { CONCAT_FILES } from './../modules/utils/main'
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
        shigella = asm.filter { cfg, contigs -> cfg.species =~ 'Shigella'}
        sonnei = reads.filter { cfg, reads -> cfg.species == 'Shigella sonnei'}
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
        // println klebs_typers.view()
        klebs_version = KLEBORATE.out.version.map {cfg, version -> version }.collect()
        ECTYPER ( ecoli )
        ecoli_typers = ECTYPER.out.typer.map {cfg, typer -> typer }.collect()
        ecoli_version = ECTYPER.out.version.map {cfg, version -> version }.collect()
        EMMTYPER ( igas )
        emm_typers = EMMTYPER.out.typer.map {cfg, typer -> typer }.collect()
        emm_version = EMMTYPER.out.version.map {cfg, version -> version }.collect()
        SHIGAPASS ( shigella )
        shigella_typers = SHIGAPASS.out.typer.map {cfg, typer -> typer }.collect()
        shigella_version = SHIGAPASS.out.version.map {cfg, version -> version }.collect()
        SONNEITYPE ( sonnei )
        sonnei_typers = SONNEITYPE.out.typer.map {cfg, typer -> typer }.collect()
        sonnei_version = SONNEITYPE.out.version.map {cfg, version -> version }.collect()

        // typers = lissero_typers.concat ( salmo_typers, nmen_typers, ngono_typers, klebs_typers, ecoli_typers,emm_typers ).flatten().toList().map { files -> tuple("typer", files)}
        // println typers.view()
        // println klebs_typers.view()
        kleb_typing = CONCAT_KLEBORATE ( klebs_typers.map {files -> tuple("kleborate", files)} )
        salmo_typing = CONCAT_SALMO ( salmo_typers.map {files -> tuple("stype", files)} )
        lissero_typing = CONCAT_LISSERO ( lissero_typers.map {files -> tuple("lissero", files)} )
        nmen_typing = CONCAT_MENINGOTYPE ( nmen_typers.map {files -> tuple("meningotype", files)} )
        ngono_typing = CONCAT_NGMASTER ( ngono_typers.map {files -> tuple("ngmaster", files)} )
        ecoli_typing = CONCAT_ECTYPER ( ecoli_typers.map {files -> tuple("ectyper", files)} )
        emm_typing = CONCAT_EMMTYPER ( emm_typers.map {files -> tuple("emmtyper", files)} )
        shigapass = CONCAT_SHIGAPASS ( shigella_typers.map {files -> tuple("shigapass", files)} )
        sonnei_typing = CONCAT_SONNEITYPE ( sonnei_typers.map {files -> tuple("sonneitype", files)} )
        collated_typers =  lissero_typing.concat(kleb_typing,salmo_typing,nmen_typing,ngono_typing,ecoli_typing,emm_typing, shigapass, sonnei_typing)
        versions = lissero_version.concat ( salmo_version, nmen_version, ngono_version, klebs_version, ecoli_version, emm_version ).map { files -> tuple("version_serotypes", files)}
        
        // CONCAT_FILES ( typers )
        CSVTK_UNIQ ( versions )
        collated_versions = CSVTK_UNIQ.out.collated
        // collated_typers = CONCAT_FILES.out.collated
    emit:
        collated_typers = collated_typers.toList()
        collated_versions

}

workflow CONCAT_KLEBORATE {
    take:
        klebs
    main:
        CONCAT_FILES ( klebs )
    emit:
        collated_klebs = CONCAT_FILES.out.collated
}

workflow CONCAT_SHIGAPASS {
    take:

        shigella
    main:
        CONCAT_FILES ( shigella )
    emit:
        collated_shigella = CONCAT_FILES.out.collated
}


workflow CONCAT_SONNEITYPE {
    take:

        sonnei
    main:
        CONCAT_FILES ( sonnei )
    emit:
        collated_sonnei = CONCAT_FILES.out.collated
}

workflow CONCAT_EMMTYPER {
    take:
        emm
    main:
        CONCAT_FILES ( emm )
    emit:
        collated_emm = CONCAT_FILES.out.collated
}


workflow CONCAT_LISSERO {
    take:
        lissero
    main:
        CONCAT_FILES ( lissero )
    emit:
        collated_lissero = CONCAT_FILES.out.collated
}


workflow CONCAT_SALMO {
    take:
        salmo
    main:
        CONCAT_FILES ( salmo )
    emit:
        collated_salmo = CONCAT_FILES.out.collated
}


workflow CONCAT_NGMASTER {
    take:
        ngmaster
    main:
        CONCAT_FILES ( ngmaster )
    emit:
        collated_ngmaster = CONCAT_FILES.out.collated
}


workflow CONCAT_MENINGOTYPE {
    take:
        meningotype
    main:
        CONCAT_FILES ( meningotype )
    emit:
        collated_meningotype = CONCAT_FILES.out.collated
}

workflow CONCAT_ECTYPER {
    take:
        ecoli
    main:
        CONCAT_FILES ( ecoli )
    emit:
        collated_ectyper = CONCAT_FILES.out.collated
}