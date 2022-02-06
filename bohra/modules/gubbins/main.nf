// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process GUBBINS {
    
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode
    
    cache 'lenient'
    conda (params.enable_conda ? 'gubbins=2.4.1 snp-sites=2.5.1 snippy=4.4.5' : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container 'https://depot.galaxyproject.org/singularity/fastp:0.20.1--h8b12597_0'
    // } else {
    //     container 'quay.io/biocontainers/fastp:0.20.1--h8b12597_0'
    // }

    input:
    path(full_core)

    output:
    path('gubbins.aln'), emit: gubbins

    // tuple val(meta), path('spades.log'), emit: log

    script:
    """
    snippy-clean_full_aln $full_core  > clean.full.aln
    run_gubbins.py -c $task.cpus  --prefix clean clean.full.aln
    snp-sites -c clean.filtered_polymorphic_sites.fasta > gubbins.aln
    """
    
}
