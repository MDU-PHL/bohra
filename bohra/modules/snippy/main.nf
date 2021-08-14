// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SNIPPY {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    
    // scratch true
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
    tuple val(meta), path('snps.aligned.fa'), emit: aln
    tuple val(meta), path('snps.raw.vcf'), emit: raw_vcf
    tuple val(meta), path('snps.vcf'), emit: vcf
    tuple val(meta), path('snps.log'), emit: log
    tuple val(meta), path('snps.tab'), emit: tab

    script:
    """
    snippy --outdir ${meta.id} --R1 ${reads[0]} --R2 ${reads[0]} --reference $launchDir/${params.reference_fasta} --force --cpus $task.cpus
    cp ${meta.id}/snps.* .
    """
    
}
