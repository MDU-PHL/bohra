// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; extract_species } from './functions'

module_dir = moduleDir + "/bin"

params.options = [:]
def options    = initOptions(params.options)


process TBTAMR {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"${meta.id}", publish_id:meta.id) }
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-tbtamr"
        } else {
            conda 'bioconda::tbtamr=1.0.3 csvtk'
        }
    } else {
        conda null
    }
    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("tbtamr_linelist_report.txt"), emit: tbtamr_txt
    tuple val(meta),path('*_stats.txt'), emit : tbtamr_stats
    tuple val(meta),path('*_variants.csv'), emit : tbtamr_variants
    tuple val(meta), path('version_tbtamr.txt'), emit: version

    script:
    
    """
    tbtamr full -1 ${reads[0]} -2 ${reads[1]} -t 4 -s ${meta.id} --call_lineage
    csvtk csv2tab ${meta.id}/tbtamr_linelist_report.csv tbtamr_linelist_report.txt
    cp ${meta.id}/*_stats.txt .
    cp ${meta.id}/*_variants.csv .
    echo -e tbtamr'\t'\$CONDA_PREFIX'\t'\$(tbtamr -v) | csvtk add-header -t -n 'tool,conda_env,version' > version_tbtamr.txt
    """
    
}
