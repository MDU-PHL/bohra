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
    
    
    
    
    if ( params.enable_conda ) {
        if (file("${params.dependency_prefix}/mash").exists()) {
            conda "${params.dependency_prefix}/mash"
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
    tuple val(meta), path('version_mash.txt'), emit: version

    script:
    """
    mash sketch -r ${reads[0]} -m 5 -k 21 -C $meta.id -o ${meta.id}
    echo -e mash'\t'\$CONDA_PREFIX'\t'\$(mash --version)'\t'${params.mash_ref} | csvtk add-header -t -n 'tool,conda_env,version,reference' > version_mash.txt
    """
    
}
