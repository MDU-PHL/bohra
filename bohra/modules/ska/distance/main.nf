// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

module_dir = moduleDir + "/bin"


process SKA_DISTANCE {

    label 'process_medium'
    publishDir "${params.outdir}/report",
        mode: params.publish_dir_mode
    
    
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/mash" : 'mash') : null) 
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-ska2"
        } else {
            conda 'ska pandas'
        }
    } else {
        conda null
    }

    cache 'lenient'
    scratch true
    
    input:
    val(merged_skf)

    output:
    path('ska_measurement.tsv'), emit: distance_long
    path('distances.tsv'), emit: matrix

    script:
    """
    ska distance -m ${params.ska_minfreq} $merged_skf > ska_measurement.tsv 
    ${module_dir}/make_matrix.py ska_measurement.tsv 
    """
        
}
