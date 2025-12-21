// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process PROKKA {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }

    cache 'lenient'
    
   
    if ( params.enable_conda ) {
        if (file("${params.dependency_prefix}/prokka").exists()) {
            conda "${params.dependency_prefix}/prokka"
        } 
    } else {
        conda null
    }

    scratch true
    
    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path('*.gff'), emit: gff
    tuple val(meta), path("${meta.id}.txt"), emit: prokka_txt
    tuple val(meta), path('version_prokka.txt'), emit: version

    script:
    """
    prokka --outdir $meta.id --prefix $meta.id --mincontiglen 500 --notrna --fast --force $contigs --cpus $task.cpus --compliant
    cp ${meta.id}/${meta.id}.gff ${meta.id}.gff
    grep -v '^##' ${meta.id}/${meta.id}.txt > ${meta.id}.txt
    echo -e prokka'\t'\$CONDA_PREFIX'\t'\$(prokka -v  2>&1)'\t'${params.prokka_ref} | csvtk add-header -t -n 'tool,conda_env,version,reference' > version_prokka.txt
    """
    
}
