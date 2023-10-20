#!/usr/bin/env nextflow
nextflow.enable.dsl=2
version = '1.0'

include { EMMTYPER } from '/../../bohra/modules/emmtyper'

process TEST_EMMTYPER {
    debug true
    input:
        tuple val(meta), val(emmtyper)
    exec:
        t = file("${params.truth}")
        myReader = emmtyper.newReader()
        result = myReader.getText()
        truthReader = t.newReader()
        truth =truthReader.getText()
        
        md5 = result.md5()
        assert md5 == truth.md5()
        
    
}

workflow {

    
    contigs = Channel.fromPath(params.contig_path)
                                                .map { file -> tuple([id: 'for_testing',single_end:false,contigs:'no_contigs'], file, 'Escherichia coli')}
    

    EMMTYPER ( contigs )

    // TEST_EMMTYPER ( EMMTYPER.out.typer )

}