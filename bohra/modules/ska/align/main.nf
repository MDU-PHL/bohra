// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SKA_ALIGN {

    label 'process_medium'
    publishDir "${params.outdir}/report",
        mode: params.publish_dir_mode
    
    
    
    
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

    output:
    path('ska.aln'), emit: aln
    
    script:
    """
    ska align -o ska.aln $merged_skf -m ${params.ska_minfreq} ${params.ska_alnargs}
    """
        
}
