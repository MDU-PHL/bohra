// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SEQTK {
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

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('seqtk_stats.txt'), emit: seqtk_stats

    script:
    """
    echo -e Bases'\t'C%'\t'G%'\t'AvgQ > seqtk_stats.txt
    cat ${reads[0]} ${reads[1]} | seqtk fqchk -q0 -  | grep ALL | cut -f2,4,5,8 >>seqtk_stats.txt
    """
    
}
