#!/usr/bin/env nextflow

include { KRAKEN2 } from './../../../modules/kraken2'


workflow {
    reads = Channel.fromFilePairs(["bohra/tests/data/*/*_{1,2}.f*q.gz"])
                .map { sample, files -> tuple([id: 'ERR1102348', single_end:false, contigs: 'no_contigs'], files)}
    
}