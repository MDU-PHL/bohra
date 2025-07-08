// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'
module_dir = moduleDir + "/bin"
params.options = [:]
def options    = initOptions(params.options)

process EXTRACT_GROUPS {
    // tag "$meta.id"
    label 'process_upper_medium'
    //  publishDir "${params.outdir}",
    //     mode: params.publish_dir_mode,
    //     saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"${meta.id}", publish_id:meta.id) }  
    cache 'lenient'
    
    
    scratch true
    // errorStrategy 'ignore'

    input:
    val(groups)
    

    output:    
    path('group.txt'),  emit: pangenome_groups

    script:
    """
    $module_dir/extract_groups.py $groups
    """
    
}

