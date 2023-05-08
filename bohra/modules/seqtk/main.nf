// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SEQTK {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    
    cache 'lenient'
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/seqtk" : 'seqtk') : null) 
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-seqtk"
        } else {
            conda 'bioconda::seqtk'
        }
    } else {
        conda null
    }
    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('seqtk_stats.txt'), emit: seqtk_stats

    script:
    """
    echo -e Bases'\t'C%'\t'G%'\t'AvgQ > seqtk_stats.txt
    cat ${reads[0]} ${reads[1]} | seqtk fqchk -q0 -  | grep ALL | cut -f2,4,5,8 >>seqtk_stats.txt
    """
    
}
