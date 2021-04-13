// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process ABRITAMR {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    
    cpus options.args2// args2 needs to be cpus for shovill
    cache 'lenient'

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path('summary_matches.txt'), emit: matches
    tuple val(meta), path('summary_partials.txt'), emit: partials
    path '*.version.txt'                  , emit: version

    script:
    // Added soft-links to original fastqs for consistent naming in MultiQC
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    abritamr -px $meta.id -c $contigs -j $task.cpus
    cp ${prefix}/summary_* .
    echo \$(abritamr --version 2>&1) | sed -e "s/abritamr //g" > ${software}.version.txt
    """
    
}
