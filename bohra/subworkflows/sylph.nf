#!/usr/bin/env nextflow
include { SYLPH } from './../modules/sylph'
include { EXTRACT_SPECIES } from './../modules/extract_species/main'



workflow RUN_SYLPH {
    take:
        sequences
    main:
        SYLPH ( sequences )
        
        EXTRACT_SPECIES ( SYLPH.out.species, "sylph" )
        // println EXTRACT_SPECIES.out.extracted_species.view()
        versions = SYLPH.out.version.map { cfg, file -> file }.collect()
                                        .map { files -> tuple("version_sylph", files) }
        CSVTK_UNIQ ( versions )
    emit:
        species_raw = SYLPH.out.species_raw
        species = SYLPH.out.species
        species_obs = EXTRACT_SPECIES.out.extracted_species
        version = CSVTK_UNIQ.out.collated
        
}