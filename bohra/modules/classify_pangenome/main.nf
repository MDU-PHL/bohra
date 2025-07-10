// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'
module_dir = moduleDir + "/bin"
params.options = [:]
def options    = initOptions(params.options)

process CLASSIFY_PANGENOME {
    // tag "$meta.id"
    label 'process_upper_medium'
    //  publishDir "${params.outdir}",
    //     mode: params.publish_dir_mode,
    //     saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"${meta.id}", publish_id:meta.id) }  
    cache 'lenient'
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-classify-pangenome"
        } else {
            conda 'conda-forge::r-data-table r-ggplot2 r-optparse'
        }
    } else {
        conda null
    }
    scratch true
    // errorStrategy 'ignore'

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

