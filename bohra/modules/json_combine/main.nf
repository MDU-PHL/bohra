// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'
module_dir = moduleDir + "/bin"
params.options = [:]
def options    = initOptions(params.options)

process JSON_COMBINE {

    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:params.report_outdir, publish_id:report_outdir) }
    

    cache 'lenient'
    scratch true
    input:
    tuple val(output_name), val(input)

    output:
    path("*.json"), emit: collated

    script:
    // Added soft-links to original fastqs for consistent naming in MultiQC
    def input_files = input.join(' ')
    """
    $module_dir/combine.py $input_files > ${output_name}.json
    """
        
}
