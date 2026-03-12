// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SEQKIT_SEQ {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    
    cache 'lenient'
    scratch true
     
    if ( params.enable_conda ) {
        if (file("${params.dependency_prefix}/seqkit").exists()) {
            conda "${params.dependency_prefix}/seqkit"
        } 
    } else {
        conda null
    }


    input:
    tuple val(meta), path(input_files)

    output:
    tuple val(meta), path('contigs.fa'), emit: contigs

    script:
    
    """
    seqkit seq -m ${params.min_contig_length} $input_files > contigs_filtered.fa
    """
       
}
