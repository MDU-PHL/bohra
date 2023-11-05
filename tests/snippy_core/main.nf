#!/usr/bin/env nextflow
nextflow.enable.dsl=2
version = '1.0'

include { SNIPPY_CORE } from './../../bohra/modules/snippy_core'

process TEST_SNIPPY_CORE {
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
    aln = Channel.fromPath("${params.aln_path}/*/*.fa").map { file -> tuple([id: file.getParent().getName(),single_end:false,contigs:'no_contigs'], file)}
    alns = aln.map { cfg, aln -> aln.getParent() }.collect() 
    reference = Channel.fromPath( "${params.reference}")

    SNIPPY_CORE ( alns,reference )  
    TEST_SNIPPY_CORE ( SNIPPY_CORE.out.core_aln )

}