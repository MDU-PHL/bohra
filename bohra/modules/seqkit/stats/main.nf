// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SEQKIT_STATS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    
    cache 'lenient'
    conda (params.enable_conda ? 'bioconda::seqkit=2.1.0' : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container 'https://depot.galaxyproject.org/singularity/fastp:0.20.1--h8b12597_0'
    // } else {
    //     container 'quay.io/biocontainers/fastp:0.20.1--h8b12597_0'
    // }

    input:
    tuple val(meta), path(input_files)

    output:
    tuple val(meta), path('*_statistics.txt'), emit: stats

    script:
    def input_type = options.args2
    if ( input_type == 'contigs'){
    """
    cat $input_files | seqkit stats --all -T  > assembly_statistics.txt
    """
    } else {   
    """
    cat ${input_files[0]} ${input_files[1]} | seqkit stats -T --all > read_statistics.txt
    """ 
    }
    
}
