#!/usr/bin/env nextflow
// include { COLLATE_KRAKEN2_ISOLATE} from './../modules/collation/main'
include { KRAKEN2 } from './../modules/kraken2'
include { EXTRACT_SPECIES } from './../modules/extract_species/main'
include { CSVTK_UNIQ } from './../modules/csvtk/main'


workflow RUN_KRAKEN {
    take:
        sequences

    main:
        
        KRAKEN2 ( sequences )
        // COLLATE_KRAKEN2_ISOLATE ( KRAKEN2.out.kraken2 )
        
        EXTRACT_SPECIES ( KRAKEN2.out.species, "kraken2" )
        versions = KRAKEN2.out.version.map { cfg, file -> file }.collect()
                                        .map { files -> tuple("version_kraken2", files) }
        CSVTK_UNIQ ( versions )
    emit:
        species_raw = KRAKEN2.out.species_raw
        species = KRAKEN2.out.species
        species_obs = EXTRACT_SPECIES.out.extracted_species
        version = CSVTK_UNIQ.out.collated
}