#!/usr/bin/env nextflow
// include { COLLATE_KRAKEN2_ISOLATE} from './../modules/collation/main'
include { KRAKEN2 } from './../modules/kraken2'
include { EXTRACT_SPECIES } from './../modules/extract_species/main'


workflow RUN_KRAKEN {
    take:
        sequences

    main:
        
        KRAKEN2 ( sequences )
        // COLLATE_KRAKEN2_ISOLATE ( KRAKEN2.out.kraken2 )
        
        EXTRACT_SPECIES ( KRAKEN2.out.species, "kraken2" )
        
    emit:
        species_raw = KRAKEN2.out.species_raw
        species = KRAKEN2.out.species
        species_obs = EXTRACT_SPECIES.out.extracted_species
}