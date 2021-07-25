#!/usr/bin/env nextflow

include { SEQTK } from './../modules/seqtk/main' 
include { COLLATE_KRAKEN2_ISOLATE; COLLATE_STATS_ISOLATE} from './../modules/collation/main'
include { SEQKIT_STATS } from './../modules/seqkit/stats/main' addParams( options: [args: '-Q -L -g', args2: 'reads'] )
include { MASH_SKETCH } from './../modules/mash/sketch/main' addParams( options: [args2: 16] )
include { MASH_TRIANGLE } from './../modules/mash/triangle/main'
include { QUICKTREE } from './../modules/quicktree/main' 
include { KRAKEN2 } from './../modules/kraken2/main' addParams( options: [args: "-db ${params.kraken2_db}", args2: 16] )

workflow READ_ANALYSIS {   

    take:
        preview
    main:
        SEQTK ( preview )
        SEQKIT_STATS ( preview )
        COMBD = SEQTK.out.seqtk_stats.join( SEQKIT_STATS.out.stats )
        COLLATE_STATS_ISOLATE ( COMBD )
        MASH_SKETCH ( preview )       
        

    emit:
        stats = COLLATE_STATS_ISOLATE.out.read_assessment
        skch = MASH_SKETCH.out.sketch

}

workflow RUN_KRAKEN {
    take:
        preview
    main:
        
        KRAKEN2 ( preview )
        COLLATE_KRAKEN2_ISOLATE ( KRAKEN2.out.kraken2 )
        
    emit:
        species = COLLATE_KRAKEN2_ISOLATE.out.species
        

}

workflow PREVIEW_NEWICK {

    take:
        sketches
        // seqkit
        // seqtk_stats
    main:
         MASH_TRIANGLE ( sketches )
         QUICKTREE ( MASH_TRIANGLE.out.mash_distances )
    emit:
        nwk = QUICKTREE.out.preveiw_tree

}

