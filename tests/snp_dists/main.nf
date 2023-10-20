#!/usr/bin/env nextflow
nextflow.enable.dsl=2
version = '1.0'

include { SNP_DISTS } from './../../modules/snp_dists'

process TEST_SNP_DISTS {
    debug true
    input:
        val dists 
    exec:
        t = file("${params.truth}")
        myReader = dists.newReader()
        result = myReader.getText()
        truthReader = t.newReader()
        truth =truthReader.getText()
        
        md5 = result.md5()
        assert md5 == truth.md5()
        
    
}

workflow {
    aln = Channel.fromPath("${params.aln_path}/core.aln")
    
    SNP_DISTS ( aln )  
    TEST_SNP_DISTS ( SNP_DISTS.out.distances )

}