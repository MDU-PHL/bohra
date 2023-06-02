#!/usr/bin/env nextflow
include { COLLATE_KRAKEN2_ISOLATE} from './../modules/collation/main'
include { KRAKEN2 } from './../modules/kraken2'
include { EXTRACT_SPECIES } from './../modules/utils/main'

workflow RUN_KRAKEN {
    take:
        preview
    main:
        
        KRAKEN2 ( preview )
        COLLATE_KRAKEN2_ISOLATE ( KRAKEN2.out.kraken2 )
        EXTRACT_SPECIES ( KRAKEN2.out.kraken2 )
        
    emit:
        kraken2 = KRAKEN2.out.kraken2
        species = COLLATE_KRAKEN2_ISOLATE.out.species
        species_obs = EXTRACT_SPECIES.out.species_obs
        

}
