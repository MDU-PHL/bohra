#!/usr/bin/env nextflow
nextflow.enable.dsl=2
version = '1.0'

include { MLST } from './../../modules/mlst'

process TEST_MLST {
    debug true
    input:
        tuple val(meta), val(mlst)
    exec:
        t = file("${params.truth}")
        myReader = mlst.newReader()
        result = myReader.getText()
        truthReader = t.newReader()
        truth =truthReader.getText()
        
        md5 = result.md5()
        assert md5 == truth.md5()
        
    
}

workflow {

    
    contigs = Channel.fromPath(params.contig_path)
                                                .map { file -> tuple([id: 'for_testing',runid:'test', barcode: '', single_end:false, job_type: '',species_exp: 'Listeria monocytogenes' ], file)}
    

    MLST ( contigs )
    println params.truth

    TEST_MLST ( MLST.out.mlst )

}