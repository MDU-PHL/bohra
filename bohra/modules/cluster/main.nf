// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'
module_dir = moduleDir + "/bin"
params.options = [:]
def options    = initOptions(params.options)

process SNP_CLUSTER {
    
    publishDir "${params.outdir}",
        mode: 'copy',
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:params.report_outdir) }
    
    if ( params.enable_conda ) {
        if (file("${params.dependency_prefix}/cluster").exists()) {
            conda "${params.dependency_prefix}/cluster"
        } else {
            conda "${moduleDir}/environment.yml"
        }
    } else {
        conda null
    }

    cache 'lenient'
    scratch true
    input:
    val(matrix) 
    
    output:
    path('clusters.txt'), emit: clusters
    

    script:

    """
    $module_dir/cluster.py $matrix ${params.cluster_threshold} ${params.cluster_method} 
    """
    
}
