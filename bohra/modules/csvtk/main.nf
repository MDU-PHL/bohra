// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process CSVTK_CONCAT {

    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    
    cache 'lenient'
    
    input:
    tuple val(output_name), val(input)

    output:
    path("*.txt"), emit: collated

    script:
    // Added soft-links to original fastqs for consistent naming in MultiQC
    def input_files = input.join(' ')
    def header = output_name =~ "mlst" ? "--no-header-row" : ""
    """
    csvtk concat $header -t -T $input_files > ${output_name}.txt
    """
        
}
