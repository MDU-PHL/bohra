// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

module_dir = moduleDir + "/bin"

params.options = [:]
def options    = initOptions(params.options)

process RUN_SMRYZR {
    
    label 'process_high'
    publishDir "${params.outdir}",
        mode: 'copy',
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"report") }
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/datasmryzr"
        } else {
            conda 'environment.yml'
        }
    } else {
        conda null
    }
    cache 'lenient'
    scratch true
    
    input:
    val(inputs)
    

    output:
    path('*.html'), emit:report
    

    script:
    def input_files = inputs.join(' ')
    """
    $module_dir/run_datasmryzr.py --job_id '${params.job_id}' \
    --mask ${params.mask} --reference ${params.reference} \
    --input_file ${launchDir}/${params.isolates} \
    --annot_cols '${params.annot_cols}' \
    --cluster_method '${params.cluster_method}' \
    --cluster_threshold '${params.cluster_threshold}' --pangenome_groups '${params.pangenome_groups}' \
    --speciation '${params.use_kraken2}' \
    --kraken2_db '${params.kraken2_db}' \
    --bkgd '${params.background_color}' \
    --text_color '${params.text_color}' \
    --results_files ${input_files} \
    """
}
