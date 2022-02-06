// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process MOBSUITE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    
    cache 'lenient'
    conda (params.enable_conda ? 'bioconda::mob_suite=3.0.3' : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container 'https://depot.galaxyproject.org/singularity/fastp:0.20.1--h8b12597_0'
    // } else {
    //     container 'quay.io/biocontainers/fastp:0.20.1--h8b12597_0'
    // }

    input:
    tuple val(meta), path(asm)

    output:
    tuple val(meta), path('contig_report.txt'), emit: contig_report
    tuple val(meta), path('mobtyper_results.txt'), emit: mobs
    tuple val(meta), path('*.fasta') optional true

    // tuple val(meta), path('spades.log'), emit: log

    script:
    """
    mob_recon -i $asm -s $meta.id -n $task.cpus -o mob
    cp mob/* .
    if [ ! -f mobtyper_results.txt ];then
        touch mobtyper_results.txt
    fi
    """
    
}
