#!/usr/bin/env nextflow
include { COLLATE_KRAKEN2_ISOLATE} from './../modules/collation/main'
include { KRAKEN2 } from './../modules/kraken2'
workflow RUN_KRAKEN {
    take:
        preview
    main:
        
        KRAKEN2 ( preview )
        COLLATE_KRAKEN2_ISOLATE ( KRAKEN2.out.kraken2 )
        
    emit:
        species = COLLATE_KRAKEN2_ISOLATE.out.species
        

}
