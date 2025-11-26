// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SNIPPY_CLEAN {
    
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode
    
    cache 'lenient'
    scratch true
    
    if ( params.enable_conda ) {
        if (file("${params.dependency_prefix}/snippy").exists()) {
            conda "${params.dependency_prefix}/snippy"
        } 
    } else {
        conda null
    }
    scratch true

    input:
    path(full_core)

    output:
    path('clean.full.aln'), emit: cleaned

    // tuple val(meta), path('spades.log'), emit: log

    script:
    """
    snippy-clean_full_aln $full_core  > clean.full.aln
    
    """
    
}
