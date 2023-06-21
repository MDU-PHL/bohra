#!/usr/bin/env nextflow
nextflow.enable.dsl=2
version = '1.0'

include { SNIPPY_CLEAN } from './../../modules/snippy_clean'

process TEST_SNIPPY_CLEAN {
    debug true
    input:
        val clean_aln 
    exec:
        t = file("${params.truth}")
        myReader = clean_aln.newReader()
        result = myReader.getText()
        truthReader = t.newReader()
        truth =truthReader.getText()
        
        md5 = result.md5()
        assert md5 == truth.md5()
        
    
}

workflow {
    aln = Channel.fromPath("${params.aln_path}/core.full.aln")
    
    SNIPPY_CLEAN ( aln )  
    TEST_SNIPPY_CLEAN ( SNIPPY_CLEAN.out.cleaned )

}