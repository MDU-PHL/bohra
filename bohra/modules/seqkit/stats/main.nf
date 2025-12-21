// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SEQKIT_STATS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    
    cache 'lenient'
    scratch true
     
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
    tuple val(meta), path('*_statistics.txt'), emit: stats
    tuple val(meta), path('version_seqkit.txt'), emit: version

    script:
    
    if ( meta.input_type == 'asm'){
    """
    cat $input_files | seqkit stats --all -T -i ${meta.id} > assembly_statistics.txt
    echo -e seqkit'\t'\$CONDA_PREFIX'\t'\$(seqkit version)'\t'${params.seqkit_ref} | csvtk add-header -t -n 'tool,conda_env,version,reference' > version_seqkit.txt
    """
    } else {   
    """
    cat ${input_files[0]} ${input_files[1]} | seqkit stats -T -i ${meta.id} --all > read_statistics.txt
    echo -e seqkit'\t'\$CONDA_PREFIX'\t'\$(seqkit version)'\t'${params.seqkit_ref} | csvtk add-header -t -n 'tool,conda_env,version,reference' > version_seqkit.txt
    """ 
    }
    
}
