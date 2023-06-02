// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process MASH_SKETCH {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    
    
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/mash" : 'mash') : null) 
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-mash"
        } else {
            conda 'mash'
        }
    } else {
        conda null
    }

    scratch true
    
    cache 'lenient'
    

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.msh'), emit: sketch

    script:
    """
    mash sketch -r ${reads[0]} -m 5 -k 21 -C $meta.id -o ${meta.id}
    """
    
}
