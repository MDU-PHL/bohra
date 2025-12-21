// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process FASTP {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode
        // saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    
    scratch true
    cache 'lenient'
    
    
    if ( params.enable_conda ) {
        if (file("${params.dependency_prefix}/fastp").exists()) {
            conda "${params.dependency_prefix}/fastp"
        } 
    } else {
        conda null
    }
    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*trim.fastq.gz'), emit: read
    tuple val(meta), path('*.json')         , emit: json
    tuple val(meta), path("version_fastp.txt"), emit: version
    

    script:
    """
    fastp \\
    --in1 ${reads[0]} \\
    --in2 ${reads[1]} \\
    --out1 ${meta.id}_1.trim.fastq.gz \\
    --out2 ${meta.id}_2.trim.fastq.gz \\
    --json ${meta.id}.fastp.json \\
    --thread $task.cpus \\
    --trim_poly_g \\
    --detect_adapter_for_pe \\
    2> ${meta.id}.fastp.log
    echo -e fastp'\t'\$CONDA_PREFIX'\t'\$(fastp -v)'\t'${params.fastp_ref} | csvtk add-header -t -n 'tool,conda_env,version,reference' > version_fastp.txt
    """
    
}
