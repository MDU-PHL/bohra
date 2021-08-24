// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process PROKKA {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }

    cache 'lenient'
    // conda (params.enable_conda ? 'bioconda::shovill=1.1.0' : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container 'https://depot.galaxyproject.org/singularity/fastp:0.20.1--h8b12597_0'
    // } else {
    //     container 'quay.io/biocontainers/fastp:0.20.1--h8b12597_0'
    // }

    scratch true
    
    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path('*.gff'), emit: gff
    tuple val(meta), path('*.txt'), emit: prokka_txt

    script:
    """
    prokka --outdir $meta.id --prefix $meta.id --mincontiglen 500 --notrna --fast --force $contigs --cpus $task.cpus
    cp ${meta.id}/${meta.id}.gff ${meta.id}.gff
    grep -v '^##' ${meta.id}/${meta.id}.txt > ${meta.id}.txt
    """
    
}