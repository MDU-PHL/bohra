// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SKA_ALIGN {

    label 'process_medium'
    publishDir "${params.outdir}/report",
        mode: params.publish_dir_mode
    
    
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/mash" : 'mash') : null) 
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-ska2"
        } else {
            conda 'mash'
        }
    } else {
        conda null
    }

    cache 'lenient'
    scratch true
    
    input:
    val(merged_skf)

    output:
    path('ska.aln'), emit: aln

    script:
    """
    ska align -o ska.aln $merged_skf -m ${params.ska_min-freq} ${params.ska_aln_args}
    """
        
}
