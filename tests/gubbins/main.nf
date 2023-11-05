#!/usr/bin/env nextflow
nextflow.enable.dsl=2
version = '1.0'

include { GUBBINS } from './../../bohra/modules/gubbins'

process TEST_GUBBINS {
    debug true
    input:
        val core_aln 
    exec:
        t = file("${params.truth}")
        myReader = core_aln.newReader()
        result = myReader.getText()
        truthReader = t.newReader()
        truth =truthReader.getText()
        
        md5 = result.md5()
        assert md5 == truth.md5()
        
    
}

workflow {
    aln = Channel.fromPath("${params.aln_path}/cleaned.full.aln")
    

    GUBBINS ( aln )  
    // TEST_GUBBINS ( GUBBINS.out.gubbins )

}