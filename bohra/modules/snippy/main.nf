// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SNIPPY {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode
        // saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    
    // scratch true
    cache 'lenient'
    conda (params.enable_conda ? 'bioconda::snippy=4.4.5' : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container 'https://depot.galaxyproject.org/singularity/fastp:0.20.1--h8b12597_0'
    // } else {
    //     container 'quay.io/biocontainers/fastp:0.20.1--h8b12597_0'
    // }

    input:
    tuple val(meta), path(reads), path(reference)

    output:
    // tuple val(meta), path("${meta.id}/*"), emit: snippy_dir
    tuple val(meta), path("${meta.id}/snps.aligned.fa"), emit: aln
    tuple val(meta), path("${meta.id}/snps.raw.vcf"), emit: raw_vcf
    tuple val(meta), path("${meta.id}/snps.vcf"), emit: vcf
    tuple val(meta), path("${meta.id}/snps.log"), emit: log
    tuple val(meta), path("${meta.id}/snps.tab"), emit: tab
    

    script:
    """
    snippy --outdir ${meta.id} \\
    --R1 ${reads[0]} --R2 ${reads[0]} \\
    --reference $reference \\
    --mapqual ${params.minmap} --basequal ${params.basequal} \\
    --mincov ${params.mincov} --minfrac ${params.minfrac} \\
    --minqual ${params.minqual} \\
    --force --cpus $task.cpus
    """
    
}
