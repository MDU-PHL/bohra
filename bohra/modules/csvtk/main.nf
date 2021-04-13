// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process CSVTK_CONCAT {

    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode
    
    cache 'lenient'
    
    input:
    tuple val(output_name), val(input)

    output:
    path("*.txt"), emit: collated

    script:
    // Added soft-links to original fastqs for consistent naming in MultiQC
    def software = getSoftwareName(task.process)
    def input_files = input.join(' ')
    """
    csvtk concat -t -T $input_files > ${output_name}.txt
    """
        
}
