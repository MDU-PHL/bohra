#!/usr/bin/env nextflow
nextflow.enable.dsl=2
version = '1.0'

include { KLEBORATE } from './../../modules/kleborate'

process TEST_KLEBORATE {
    debug true
    input:
        tuple val(meta), val(kleborate)
    exec:
        t = file("${params.truth}")
        myReader = kleborate.newReader()
        result = myReader.getText()
        truthReader = t.newReader()
        truth =truthReader.getText()
        
        md5 = result.md5()
        assert md5 == truth.md5()
        
    
}

workflow {

    
    contigs = Channel.fromPath(params.contig_path)
                                                .map { file -> tuple([id: 'for_testing',single_end:false,contigs:'no_contigs'], file, 'Klebsiella oxytoca')}
    

    KLEBORATE ( contigs )

    // TEST_KLEBORATE ( KLEBORATE.out.typer )

}