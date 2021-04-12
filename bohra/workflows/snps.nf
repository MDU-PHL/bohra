#!/usr/bin/env nextflow

include { SNIPPY } from './../modules/snippy/main' 
// include { COLLATE_KRAKEN2_ISOLATE; COLLATE_STATS_ISOLATE} from './../modules/collation/main'
// include { SEQKIT_STATS } from './../modules/seqkit/stats/main' addParams( options: [args: '-Q -L -g', args2: 'reads'] )
// include { MASH_SKETCH } from './../modules/mash/sketch/main' addParams( options: [args2: 16] )
// include { MASH_TRIANGLE } from './../modules/mash/triangle/main'
// include { QUICKTREE } from './../modules/quicktree/main' 
// include { KRAKEN2 } from './../modules/kraken2/main' addParams( options: [args: "-db ${params.kraken2_db}", args2: 16] )

workflow RUN_SNIPPY {   

    take:
        reads
    main:
        SNIPPY ( reads )       

}

