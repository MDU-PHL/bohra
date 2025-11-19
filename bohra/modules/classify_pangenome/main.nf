// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'
module_dir = moduleDir + "/bin"
params.options = [:]
def options    = initOptions(params.options)

process CLASSIFY_PANGENOME {
    
    label 'process_medium'
    cache 'lenient'
    
    if ( params.enable_conda ) {
        if (file("${params.dependency_prefix}/classify-pangenome").exists()) {
            conda "${params.dependency_prefix}/classify-pangenome"
        } else {
            conda "${moduleDir}/environment.yml"
        }
    } else {
        conda null
    }
    scratch true
    errorStrategy 'ignore'
    
    
    input:
    path(rtab)
    path(groups)

    output:    
    path('classification.tab'),  emit: pangenome_classification

    script:
    """
    $module_dir/classify_genes.R \\
        -p $rtab \\
        -g $groups \\
        -s 1 
    cp out/classification.tab classification.tab
    """
    
}

