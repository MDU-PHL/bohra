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
        
    emit:
        species_raw = SYLPH.out.species_raw
        species = SYLPH.out.species
        species_obs = EXTRACT_SPECIES.out.extracted_species
        
}