// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process MLST {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    
    cache 'lenient'
    
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/mlst" : 'bioconda::mlst=2.19.0') : null) 
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-mlst"
        } else {
            conda 'mlst=2.19.0'
        }
    } else {
        conda null
    }

    scratch true

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path('mlst.txt'), emit: mlst
    tuple val(meta), path('mlst.json'), emit: json

    script:
    def _blast_db = params.blast_db != 'no_db' ? "--blastdb ${params.blast_db}" : ""
    def _publst_db = params.data_dir != 'no_db' ? "--datadir ${params.data_dir}" : ""
    def exclude = params.mlst_exclude != '' ? "--exclude ${params.mlst_exclude}" : ""
    """
    mlst --json mlst.json --label $meta.id --nopath $contigs $_blast_db $_publst_db  $exclude > mlst.txt
    """
    
}
