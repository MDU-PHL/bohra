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
    // scratch true
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/seqkit" : 'csvtk seqkit=2.1.0') : null) 
    if ( params.enable_conda ) {
        if (file("${params.dependency_prefix}/seqkit").exists()) {
            conda "${params.dependency_prefix}/seqkit"
        } 
    } else {
        conda null
    }
    

    input:
    tuple val(meta), path(input_files)

    output:
    tuple val(meta), path('read_qual.txt'), emit: stats
    tuple val(meta), path('version_seqkit.txt'), emit: version
    script:
    """
    cat ${input_files[0]} ${input_files[1]} \
    | seqkit fx2tab -H --name --only-id --avg-qual --gc \
    | csvtk summary -t -i -f 2:mean,3:mean > read_qual.txt
    echo -e seqkit'\t'\$CONDA_PREFIX'\t'\$(seqkit version)'\t'${params.seqkit_ref} | csvtk add-header -t -n 'tool,conda_env,version,reference' > version_seqkit.txt
    """

    
    
}

