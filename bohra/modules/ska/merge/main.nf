// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SKA_MERGE {
    
    label 'process_medium'
    publishDir "${params.outdir}/report",
        mode: params.publish_dir_mode
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}/bohra-ska2").exists()) {
            conda "${params.conda_path}/bohra-ska2"
        } else {
            conda "${moduleDir}/environment.yml"
        }
    } else {
        conda null
    }

    scratch true
    
    cache 'lenient'
    

    input:
    val(skf)

    output:
    path("merged.skf"), emit: merged_skf

    script:
    def skfs = skf.join(' ')
    """
    ska merge -o merged.skf $skfs
    """
    
}
