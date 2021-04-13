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
    

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path('mlst.tab'), emit: mlst
    tuple val(meta), path('mlst.json'), emit: json
    path '*.version.txt'                  , emit: version

    script:
    // Added soft-links to original fastqs for consistent naming in MultiQC
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    mlst --json mlst.json --nopath $contigs > mlst.tab
    echo \$(mlst --version 2>&1) | sed -e "s/mlst //g" > ${software}.version.txt
    """
    
}
