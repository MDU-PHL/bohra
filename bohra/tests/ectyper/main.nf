#!/usr/bin/env nextflow
nextflow.enable.dsl=2
version = '1.0'

include { ECTYPER } from './../../modules/ectyper'

process TEST_ECTYPER {
    debug true
    input:
        tuple val(meta), val(ectyper)
    exec:
        t = file("${params.truth}")
        myReader = ectyper.newReader()
        result = myReader.getText()
        truthReader = t.newReader()
        truth =truthReader.getText()
        
        md5 = result.md5()
        assert md5 == truth.md5()
        
    
}

workflow {

    
    contigs = Channel.fromPath(params.contig_path)
                                                .map { file -> tuple([id: 'for_testing',single_end:false,contigs:'no_contigs'], file, 'Escherichia coli')}
    

    ECTYPER ( contigs )

    // TEST_ECTYPER ( ECTYPER.out.typer )

}