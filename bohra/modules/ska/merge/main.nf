// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SKA_MERGE {
    
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
