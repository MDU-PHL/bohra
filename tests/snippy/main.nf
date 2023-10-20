#!/usr/bin/env nextflow
nextflow.enable.dsl=2
version = '1.0'

include { SNIPPY } from '/../../bohra/modules/snippy'

process TEST_SNIPPY {
    debug true
    input:
        tuple val(meta), val(tab)
    exec:
        t = file("${params.truth}")
        myReader = tab.newReader()
        result = myReader.getText()
        truthReader = t.newReader()
        truth =truthReader.getText()
        
        md5 = result.md5()
        assert md5 == truth.md5()
        
    
}

workflow {

    
    reads = Channel.fromFilePairs(["${params.read_path}/*/R{1,2}*.f*q.gz","${params.read_path}/*/*_{1,2}.f*q.gz"])
                .map { sample, files -> tuple([id: 'for_testing', single_end:false, contigs: 'no_contigs'], files)}

    reference = Channel.fromPath( "${params.reference}")

    SNIPPY ( reads.combine( reference ) )  

    TEST_SNIPPY ( SNIPPY.out.tab )

}