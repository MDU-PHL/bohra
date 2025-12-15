// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

module_dir = moduleDir + "/bin"

params.options = [:]
def options    = initOptions(params.options)

process CHECK_CORE {
    // tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: 'copy',
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:params.report_outdir) }
    

    cache 'lenient'
    scratch true
    
    input:
    val(snippy_qc)

    output:
    path("core_genome_stats.txt"), emit: stats
    
    // def cre = snippy_qc.join(' ')
    script:
    """
    $module_dir/check_core.py --snippy_qc ${snippy_qc} --output core_genome_stats.txt 
    """
    
}


process FILTER_CORE {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: 'copy',
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:params.report_outdir) }
    

    cache 'lenient'
    scratch true
    
    input:
    tuple val(meta), val(aln), path(core_qc)

    output:
    
    tuple val(meta), stdout, emit: aln_filter
    
    // def core = snippy_qc.join(' ')
    // for now this will add exclude to every aln if any are outliers - this behviour may change in the future
    script:
    """
    $module_dir/filter.py --qc ${core_qc} --strict ${params.ignore_warnings} --seqid ${meta.id}
    """
    
}
