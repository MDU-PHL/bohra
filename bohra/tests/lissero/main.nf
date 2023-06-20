#!/usr/bin/env nextflow
nextflow.enable.dsl=2
version = '1.0'

include { LISSERO } from './../../modules/lissero'

process TEST_LISSERO {
    debug true
    input:
        tuple val(meta), val(lissero)
    exec:
        t = file("${params.truth}")
        myReader = lissero.newReader()
        result = myReader.getText()
        truthReader = t.newReader()
        truth =truthReader.getText()
        
        md5 = result.md5()
        assert md5 == truth.md5()
        
    
}

workflow {

    
    contigs = Channel.fromPath(params.contig_path)
                                                .map { file -> tuple([id: 'for_testing',single_end:false,contigs:'no_contigs'], file, 'Listeria monocytogenes')}
    

    LISSERO ( contigs )

    TEST_LISSERO ( LISSERO.out.typer )

}