// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'
module_dir = moduleDir + "/bin"
params.options = [:]
def options    = initOptions(params.options)

process CLASSIFY_PANGENOME {
    
    label 'process_upper_medium'
    cache 'lenient'
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}/bohra-classify-pangenome").exists()) {
            conda "${params.conda_path}/bohra-classify-pangenome"
        } else {
            conda "${moduleDir}/environment.yml"
        }
    } else {
        conda null
    }
    scratch true
    
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

