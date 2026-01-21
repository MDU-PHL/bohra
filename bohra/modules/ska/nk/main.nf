// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

module_dir = moduleDir + "/bin"


process SKA_NK {

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

    output:
    path('ska_summary_stats.tsv'), emit: summary_stats
    


    script:
    """
    ska nk $merged_skf > tmp_stats.txt
    ${module_dir}/summarise.py tmp_stats.txt ska_summary_stats.tsv
    """
        
}
