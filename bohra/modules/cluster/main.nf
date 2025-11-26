// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'
module_dir = moduleDir + "/bin"
params.options = [:]
def options    = initOptions(params.options)

process SNP_CLUSTER {
    
    label 'process_medium'
    publishDir "${params.outdir}/report",
        mode: params.publish_dir_mode
    
    // if ( params.enable_conda ) {
    //     if (file("${params.dependency_prefix}/cluster").exists()) {
    //         conda "${params.dependency_prefix}/cluster"
    //     } else {
    //         conda "${moduleDir}/environment.yml"
    //     }
    // } else {
    //     conda null
    // }

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
