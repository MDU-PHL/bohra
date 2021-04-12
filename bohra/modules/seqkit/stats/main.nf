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
    // conda (params.enable_conda ? 'bioconda::shovill=1.1.0' : null)
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
    // Added soft-links to original fastqs for consistent naming in MultiQC
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def input_type = options.args2
    if ( input_type == 'contigs'){
    """
    seqkit stats -T $input_files > assembly_statistics.txt
    """
    } else {   
    """
    [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${input_files[0]} ${prefix}_1.fastq.gz
    [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${input_files[1]} ${prefix}_2.fastq.gz
    seqkit stats -T ${prefix}_1.fastq.gz  ${prefix}_2.fastq.gz > read_statistics.txt
    """ 
    }
    
}
