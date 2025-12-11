// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

module_dir = moduleDir + "/bin"

params.options = [:]
def options    = initOptions(params.options)

process CHECK_CORE {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: 'copy',
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:params.report_outdir) }
    

    cache 'lenient'
    scratch true
    
    input:
    path(snippy_qc)

    output:
    
    path("core_genome_stats.txt"), emit: stats
    
    def core = snippy_qc.join(' ')
    script:
    """
    $module_dir/check_core.py --snippy_qc ${core} --output core_genome_stats.txt
    """
    
}
