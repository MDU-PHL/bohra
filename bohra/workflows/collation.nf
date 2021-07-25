#!/usr/bin/env nextflow


include { COLLATE_KRAKENS; COLLATE_SEQDATA; COMPILE } from './../modules/collation/main'


workflow COLLATE_KRAKEN {
    take:
        species
    main:
        COLLATE_KRAKENS ( species )
        
    emit:
        collated_species = COLLATE_KRAKENS.out.collated_species
        

}

workflow COLLATE_SEQS {
    take:
        stats
    main:
        COLLATE_SEQDATA ( stats )
        
    emit:
        collated_seqdata = COLLATE_SEQDATA.out.collated_seqdata
        

}


workflow WRITE_HTML {
    take:
        results
    main:
        COMPILE ( results )
        
    emit:
        written_html = COMPILE.out.html
        

}