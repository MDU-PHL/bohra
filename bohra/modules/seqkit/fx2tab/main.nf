// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SEQKIT_GC {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    
    cache 'lenient'
    
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/seqkit" : 'csvtk seqkit=2.1.0') : null) 
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-seqkit"
        } else {
            conda 'csvtk seqkit=2.1.0'
        }
    } else {
        conda null
    }
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container 'https://depot.galaxyproject.org/singularity/fastp:0.20.1--h8b12597_0'
    // } else {
    //     container 'quay.io/biocontainers/fastp:0.20.1--h8b12597_0'
    // }

    input:
    tuple val(meta), path(input_files)

    output:
    tuple val(meta), path('read_qual.txt'), emit: stats

    script:
    """
    cat ${input_files[0]} ${input_files[1]} \
    | seqkit fx2tab -H --name --only-id --avg-qual --gc \
    | csvtk summary -t -i -f 2:mean,3:mean > read_qual.txt
    """

    
    
}
