#!/usr/bin/env nextflow

include {CSVTK_CONCAT } from './../modules/csvtk/main'
include { VERSION_MASH;VERSION_QUICKTREE;VERSION_PROKKA;VERSION_PANAROO;VERSION_NGMASTER;VERSION_MOBSUITE;VERSION_MENINGOTYPE;VERSION_LISSERO;VERSION_KLEBORATE;VERSION_EMMTYPER;VERSION_ECTYPER;VERSION_ANY2FASTA;VERSION_MLST;VERSION_ABRITAMR;VERSION_SHOVILL;VERSION_SPADES;VERSION_SKESA;VERSION_SNIPPY;VERSION_SNPDISTS;VERSION_SEQKIT;VERSION_KMC;VERSION_KRAKEN2;VERSION_GUBBINS; VERSION_IQTREE} from './../modules/utils/main'


workflow PREVIEW_VERSIONS {
    
    main:
        VERSION_QUICKTREE()
        quicktree = VERSION_QUICKTREE.out.version
        VERSION_MASH()
        mash = VERSION_MASH.out.version
        VERSION_SEQKIT()
        seqkit = VERSION_SEQKIT.out.version
        VERSION_KMC()
        kmc = VERSION_KMC.out.version
        versions = quicktree.concat ( mash, seqkit, kmc )
        if ( params.run_kraken ){
            VERSION_KRAKEN2()
            kraken2 = VERSION_KRAKEN2.out.version
            versions = versions.concat(kraken2)
        }
        CSVTK_CONCAT ( versions.collect().map { files -> tuple("versions", files)})
    emit:
        versions = CSVTK_CONCAT.out.collated
}


workflow SNPS_VERSIONS {
    
    main:
        VERSION_SNIPPY()
        snippy = VERSION_SNIPPY.out.version
        VERSION_SNPDISTS()
        VERSION_SEQKIT()
        VERSION_KMC()
        snp_dists = VERSION_SNPDISTS.out.version
        seqkit = VERSION_SEQKIT.out.version
        kmc = VERSION_KMC.out.version
        versions = snippy.concat ( snp_dists, seqkit, kmc )
        if ( params.run_kraken ){
            VERSION_KRAKEN2()
            kraken2 = VERSION_KRAKEN2.out.version
            versions = versions.concat(kraken2)
        }
        if ( params.gubbins ) {
            VERSION_GUBBINS()
            gubbins = VERSION_GUBBINS.out.version
            versions = versions.concat ( gubbins )
        } 
        if (params.run_iqtree ){
            VERSION_IQTREE()
            iqtree = VERSION_IQTREE.out.version
            versions = versions.concat(iqtree)
        }

        CSVTK_CONCAT ( versions.collect().map { files -> tuple("versions", files)})
    emit:
        versions = CSVTK_CONCAT.out.collated
}

workflow AMR_TYPING_VERSIONS {
    
    main:
        VERSION_SEQKIT()
        VERSION_KMC()
        seqkit = VERSION_SEQKIT.out.version
        kmc = VERSION_KMC.out.version
        versions = seqkit.concat( kmc )
        if ( params.assembler == 'shovill'){
        VERSION_SHOVILL()
        asm = VERSION_SHOVILL.out.version    
        } 
        else if ( params.assembler == 'spades' ){
        VERSION_SPADES()
        asm = VERSION_SPADES.out.version    
        } else if (params.assembler == 'skesa' ) {
        VERSION_SKESA()
        asm = VERSION_SKESA.out.version    
        }
        versions = versions.concat ( asm )

        if ( params.run_kraken ){
            VERSION_KRAKEN2()
            kraken2 = VERSION_KRAKEN2.out.version
            versions = versions.concat(kraken2)
        }
        VERSION_ABRITAMR()
        abritamr = VERSION_ABRITAMR.out.version
        VERSION_MLST()
        mlst = VERSION_MLST.out.version
        VERSION_ECTYPER()
        ectyper = VERSION_ECTYPER.out.version
        VERSION_EMMTYPER()
        emmtyper = VERSION_EMMTYPER.out.version
        VERSION_KLEBORATE()
        kleborate = VERSION_KLEBORATE.out.version
        VERSION_LISSERO()
        lissero = VERSION_LISSERO.out.version
        VERSION_MENINGOTYPE()
        meningotype = VERSION_MENINGOTYPE.out.version
        VERSION_MOBSUITE()
        mobsuite = VERSION_MOBSUITE.out.version
        VERSION_NGMASTER()
        ngmaster = VERSION_NGMASTER.out.version
        versions = versions.concat ( ngmaster, mobsuite, lissero, kleborate, emmtyper, mlst, abritamr )
        CSVTK_CONCAT ( versions.collect().map { files -> tuple("versions", files)})
    emit:
        versions = CSVTK_CONCAT.out.collated
}

workflow ASSEMBLE_VERSIONS {

    main:
        VERSION_SEQKIT()
        VERSION_KMC()
        seqkit = VERSION_SEQKIT.out.version
        kmc = VERSION_KMC.out.version
        versions = seqkit.concat( kmc )
        if ( params.assembler == 'shovill'){
        VERSION_SHOVILL()
        asm = VERSION_SHOVILL.out.version    
        } 
        else if ( params.assembler == 'spades' ){
        VERSION_SPADES()
        asm = VERSION_SPADES.out.version    
        } else if (params.assembler == 'skesa' ) {
        VERSION_SKESA()
        asm = VERSION_SKESA.out.version    
        }
        versions = versions.concat ( asm )

        if ( params.run_kraken ){
            VERSION_KRAKEN2()
            kraken2 = VERSION_KRAKEN2.out.version
            versions = versions.concat(kraken2)
        }
        VERSION_MOBSUITE()
        mobsuite = VERSION_MOBSUITE.out.version
        VERSION_PROKKA()
        prokka = VERSION_PROKKA.out.version
        versions = versions.concat( mobsuite,prokka )
        CSVTK_CONCAT ( versions.collect().map { files -> tuple("versions", files)})
    emit:
        versions = CSVTK_CONCAT.out.collated
}

workflow DEFAULT_VERSIONS {
    
    main:
        VERSION_SNIPPY()
        snippy = VERSION_SNIPPY.out.version
        VERSION_SNPDISTS()
        VERSION_SEQKIT()
        VERSION_KMC()
        snp_dists = VERSION_SNPDISTS.out.version
        seqkit = VERSION_SEQKIT.out.version
        kmc = VERSION_KMC.out.version
        if ( params.assembler == 'shovill'){
        VERSION_SHOVILL()
        asm = VERSION_SHOVILL.out.version    
        } 
        else if ( params.assembler == 'spades' ){
        VERSION_SPADES()
        asm = VERSION_SPADES.out.version    
        } else if (params.assembler == 'skesa' ) {
        VERSION_SKESA()
        asm = VERSION_SKESA.out.version    
        }
        versions = snippy.concat ( snp_dists, seqkit, kmc, asm )

        if ( params.run_kraken ){
            VERSION_KRAKEN2()
            kraken2 = VERSION_KRAKEN2.out.version
            versions = versions.concat(kraken2)
        }
        if ( params.gubbins ) {
            VERSION_GUBBINS()
            gubbins = VERSION_GUBBINS.out.version
            versions = versions.concat ( gubbins )
        } 
        if (params.run_iqtree ){
            VERSION_IQTREE()
            iqtree = VERSION_IQTREE.out.version
            versions = versions.concat(iqtree)
        }
        VERSION_ABRITAMR()
        abritamr = VERSION_ABRITAMR.out.version
        VERSION_MLST()
        mlst = VERSION_MLST.out.version
        VERSION_ANY2FASTA()
        any2fasta = VERSION_ANY2FASTA.out.version
        VERSION_ECTYPER()
        ectyper = VERSION_ECTYPER.out.version
        VERSION_EMMTYPER()
        emmtyper = VERSION_EMMTYPER.out.version
        VERSION_KLEBORATE()
        kleborate = VERSION_KLEBORATE.out.version
        VERSION_LISSERO()
        lissero = VERSION_LISSERO.out.version
        VERSION_MENINGOTYPE()
        meningotype = VERSION_MENINGOTYPE.out.version
        VERSION_MOBSUITE()
        mobsuite = VERSION_MOBSUITE.out.version
        VERSION_NGMASTER()
        ngmaster = VERSION_NGMASTER.out.version
        VERSION_PROKKA()
        prokka = VERSION_PROKKA.out.version
        versions = versions.concat ( prokka, ngmaster, mobsuite, lissero, kleborate, emmtyper, mlst, abritamr )
        CSVTK_CONCAT ( versions.collect().map { files -> tuple("versions", files)})
    emit:
        versions = CSVTK_CONCAT.out.collated
}

workflow FULL_VERSIONS {
    
    main:
        VERSION_SNIPPY()
        snippy = VERSION_SNIPPY.out.version
        VERSION_SNPDISTS()
        VERSION_SEQKIT()
        VERSION_KMC()
        snp_dists = VERSION_SNPDISTS.out.version
        seqkit = VERSION_SEQKIT.out.version
        kmc = VERSION_KMC.out.version
        if ( params.assembler == 'shovill'){
        VERSION_SHOVILL()
        asm = VERSION_SHOVILL.out.version    
        } 
        else if ( params.assembler == 'spades' ){
        VERSION_SPADES()
        asm = VERSION_SPADES.out.version    
        } else if (params.assembler == 'skesa' ) {
        VERSION_SKESA()
        asm = VERSION_SKESA.out.version    
        }
        versions = snippy.concat ( snp_dists, seqkit, kmc, asm )

        if ( params.run_kraken ){
            VERSION_KRAKEN2()
            kraken2 = VERSION_KRAKEN2.out.version
            versions = versions.concat(kraken2)
        }
        if ( params.gubbins ) {
            VERSION_GUBBINS()
            gubbins = VERSION_GUBBINS.out.version
            versions = versions.concat ( gubbins )
        } 
        if (params.run_iqtree ){
            VERSION_IQTREE()
            iqtree = VERSION_IQTREE.out.version
            versions = versions.concat(iqtree)
        }
        VERSION_ABRITAMR()
        abritamr = VERSION_ABRITAMR.out.version
        VERSION_MLST()
        mlst = VERSION_MLST.out.version
        VERSION_ANY2FASTA()
        any2fasta = VERSION_ANY2FASTA.out.version
        VERSION_ECTYPER()
        ectyper = VERSION_ECTYPER.out.version
        VERSION_EMMTYPER()
        emmtyper = VERSION_EMMTYPER.out.version
        VERSION_KLEBORATE()
        kleborate = VERSION_KLEBORATE.out.version
        VERSION_LISSERO()
        lissero = VERSION_LISSERO.out.version
        VERSION_MENINGOTYPE()
        meningotype = VERSION_MENINGOTYPE.out.version
        VERSION_PROKKA()
        prokka = VERSION_PROKKA.out.version
        VERSION_MOBSUITE()
        mobsuite = VERSION_MOBSUITE.out.version
        VERSION_NGMASTER()
        ngmaster = VERSION_NGMASTER.out.version
        VERSION_PANAROO()
        panaroo = VERSION_PANAROO.out.version
        versions = versions.concat ( panaroo,prokka,ngmaster, mobsuite, lissero, kleborate, emmtyper, mlst, abritamr )
        CSVTK_CONCAT ( versions.collect().map { files -> tuple("versions", files)})
    emit:
        versions = CSVTK_CONCAT.out.collated
}

