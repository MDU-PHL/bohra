// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

module_dir = moduleDir + "/bin"


process SKA_DISTANCE {

    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:params.report_outdir, publish_id:params.report_outdir) }
    
    
    
    if ( params.enable_conda ) {
        if (file("${params.dependency_prefix}/ska").exists()) {
            conda "${params.dependency_prefix}/ska"
        } 
    } else {
        conda null
    }

    cache 'lenient'
    scratch true
    
    input:
    val(merged_skf)
    val(ska_nk)

    output:
    path('ska_distance_stats.tsv'), emit: distance_long
    path('distances.tsv'), emit: matrix


    script:
    """
    ska distance -m ${params.ska_minfreq} $merged_skf > ska_distance_stats.tsv 
    ${module_dir}/make_matrix.py ska_distance_stats.tsv $ska_nk
    """
        
}
