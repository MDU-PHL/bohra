#!/usr/bin/env nextflow

include { PANAROO } from './../modules/panaroo/main' 


workflow RUN_PANAROO {   

    take:
        gff
    main:
        PANAROO ( gff )     

        
    emit:
        
        roary = PANAROO.out.pangenome_summary
        version = PANAROO.out.version

}