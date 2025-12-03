// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'
module_dir = moduleDir + "/bin"
params.options = [:]
def options    = initOptions(params.options)

process CONCAT_FILES {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:params.report_outdir, publish_id:report_outdir) }
    
    errorStrategy 'ignore'
    scratch true

    cache 'lenient'
    
    input:
    tuple val(output_name), val(input)

    output:
    path("*.txt"), emit: collated

    script:
    // Added soft-links to original fastqs for consistent naming in MultiQC
    def input_files = input.join(' ')
    """
    $module_dir/concat.py ${output_name} $input_files 
    """ 
        
}


process BOHRA_VERSION {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:params.report_outdir, publish_id:params.report_outdir) }
    
    errorStrategy 'ignore'
    scratch true

    cache 'lenient'

    output:
    path("version_bohra.txt"), emit: collated

    script:
    
    """
    echo -e bohra'\t'\$CONDA_PREFIX'\t'\$(bohra --version | cut -f 3 -d' ') | csvtk space2tab | csvtk add-header -t -n 'tool,conda_env,version' > version_bohra.txt
    """ 
        
}