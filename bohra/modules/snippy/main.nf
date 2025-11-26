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
    
    scratch true
    cache 'lenient'
    
    
    if ( params.enable_conda ) {
        if (file("${params.dependency_prefix}/snippy").exists()) {
            conda "${params.dependency_prefix}/snippy"
        } 
    } else {
        conda null
    }
    input:
    tuple val(meta), path(reads), path(reference)

    output:
    tuple val(meta), path("${meta.id}/snps.aligned.fa"), emit: aln
    tuple val(meta), path("${meta.id}/snps.raw.vcf"), emit: raw_vcf
    tuple val(meta), path("${meta.id}/snps.vcf"), emit: vcf
    tuple val(meta), path("${meta.id}/snps.log"), emit: log
    tuple val(meta), path("${meta.id}/snps.tab"), emit: tab
    tuple val(meta), path("version_snippy.txt"), emit: version
    

    script:
    """
    snippy --outdir ${meta.id} \\
    --R1 ${reads[0]} --R2 ${reads[1]} \\
    --reference $reference \\
    --mapqual ${params.minmap} --basequal ${params.basequal} \\
    --mincov ${params.mincov} --minfrac ${params.minfrac} \\
    --minqual ${params.minqual} \\
    --force --cpus $task.cpus
    echo -e snippy'\t'\$CONDA_PREFIX'\t'\$(snippy -v 2>&1)'\t'${params.snippy_ref} | csvtk add-header -t -n 'tool,conda_env,version,reference' > version_snippy.txt
    """
    
}
