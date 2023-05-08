// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SNIPPY_CLEAN {
    
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode
    
    cache 'lenient'
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/snippy" : 'snippy=4.4.5') : null)
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-snippy"
        } else {
            conda 'bioconda::snippy=4.4.5'
        }
    } else {
        conda null
    }


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
