// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'
module_dir = moduleDir + "/bin"
params.options = [:]
def options    = initOptions(params.options)

process SNP_CLUSTER {
    
    label 'process_medium'
    publishDir "${params.outdir}/report",
        mode: params.publish_dir_mode
    
    cache 'lenient'
    
    input:
    val(matrix) // this needs to be a list of sample! not .aln since snippy core uses relative path and the name of the folder to name results!
    
    output:
    path('clusters_*.txt'), emit: clusters
    

    script:

    """
    $module_dir/cluster.py $matrix ${params.snpthresholds} ${params.linkage} ${params.jobname}
    """
    
}
