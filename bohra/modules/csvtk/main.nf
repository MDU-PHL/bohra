// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process CSVTK_CONCAT {

    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/csvtk" : 'bioconda::csvtk') : null)
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-csvtk"
        } else {
            conda 'bioconda::csvtk'
        }
    } else {
        conda null
    }

    cache 'lenient'
    
    input:
    tuple val(output_name), val(input)

    output:
    path("*.txt"), emit: collated

    script:
    // Added soft-links to original fastqs for consistent naming in MultiQC
    def input_files = input.join(' ')
    """
    csvtk concat -t -T $input_files > ${output_name}.txt
    """
        
}
