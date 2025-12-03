// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process CSVTK_CONCAT {

    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:params.report_outdir, publish_id:report_outdir) }
    
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/csvtk" : 'bioconda::csvtk') : null)
    
    // if ( params.enable_conda ) {
    //     if (file("${params.dependency_prefix}/csvtk").exists()) {
    //         conda "${params.dependency_prefix}/csvtk"
    //     } else {
    //         conda "${moduleDir}/environment.yml"
    //     }
    // } else {
    //     conda null
    // }

    cache 'lenient'
    scratch true
    input:
    tuple val(output_name), val(input)

    output:
    path("*.txt"), emit: collated

    script:
    // Added soft-links to original fastqs for consistent naming in MultiQC
    def input_files = input.join(' ')
    """
    csvtk -t concat -k -u '' $input_files > ${output_name}.txt
    """
        
}


process CSVTK_UNIQ {

    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:params.report_outdir, publish_id:report_outdir) }
    
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/csvtk" : 'bioconda::csvtk') : null)
    
    // if ( params.enable_conda ) {
    //     if (file("${params.dependency_prefix}/csvtk").exists()) {
    //         conda "${params.dependency_prefix}/csvtk"
    //     } else {
    //         conda "${moduleDir}/environment.yml"
    //     }
    // } else {
    //     conda null
    // }

    cache 'lenient'
    
    input:
    tuple val(output_name), val(input)

    output:
    path("*.txt"), emit: collated

    script:
    // Added soft-links to original fastqs for consistent naming in MultiQC
    def input_files = input.join(' ')
    """
    csvtk concat -u '' -t $input_files | csvtk grep -t -v -f 2 -r -p 'Not Applicable' | csvtk uniq -t > ${output_name}.txt
    """
        
}
