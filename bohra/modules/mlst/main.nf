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
    // conda (params.enable_conda ? 'bioconda::mlst=2.19.0' : null)
    scratch true

    input:
    tuple val(meta), path(contigs), path(blast_db), path(pubmlst_db)

    output:
    tuple val(meta), path('mlst.csv'), emit: mlst
    tuple val(meta), path('mlst.json'), emit: json

    script:
    // def _blast_db = blast_db  ? "--blastdb ${blast_db}" : ""
    // def _datadir = data_dir  ? "--datadir ${data_dir}" : ""
    """
    mlst --csv --json mlst.json --label $meta.id --nopath $contigs --blastdb $params.blast_db --datadir $params.data_dir > mlst.csv
    """
    
}
