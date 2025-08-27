// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'
module_dir = moduleDir + "/bin"
params.options = [:]
def options    = initOptions(params.options)

process EXTRACT_GROUPS {
    // tag "$meta.id"
    label 'process_medium'
    cache 'lenient'
    
    
    scratch true
    
    input:
    val(groups)
    

    output:    
    path('group.txt'),  emit: pangenome_groups

    script:
    """
    $module_dir/extract_groups.py $groups
    """
    
}

