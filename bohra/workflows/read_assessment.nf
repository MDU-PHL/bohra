

include { SEQTK } from './../modules/seqtk/main' 
include { COLLATE_STATS_ISOLATE} from './../modules/collation/main'
include { SEQKIT_STATS } from './../modules/seqkit/stats/main' addParams( options: [args: '-Q -L -g', args2: 'reads'] )
include { SEQKIT_GC } from './../modules/seqkit/fx2tab/main' 
include { KMC } from './../modules/kmc/main' 
include { MASH_SKETCH } from './../modules/mash/sketch/main' 
include { MASH_TRIANGLE } from './../modules/mash/triangle/main'
include { QUICKTREE } from './../modules/quicktree/main' 

workflow READ_ANALYSIS {   

    take:
        preview
        
    main:
        
        SEQKIT_STATS ( preview )
        SEQKIT_GC ( preview )
        KMC ( preview )
        COMBD = SEQKIT_STATS.out.stats.join( SEQKIT_GC.out.stats )
        COMBD = COMBD.join( KMC.out.genome_size )
        COLLATE_STATS_ISOLATE ( COMBD )
    emit:
        stats = COLLATE_STATS_ISOLATE.out.read_assessment
        

}
