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
    
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/snippy" : 'bioconda::snippy=4.4.5') : null) 
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-snippy"
        } else {
            conda 'bioconda::snippy=4.4.5'
        }
    } else {
        conda null
    }
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
    --R1 ${reads[0]} --R2 ${reads[1]} \\
    --reference $reference \\
    --mapqual ${params.minmap} --basequal ${params.basequal} \\
    --mincov ${params.mincov} --minfrac ${params.minfrac} \\
    --minqual ${params.minqual} \\
    --force --cpus $task.cpus
    """
    
}
