// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

module_dir = moduleDir + "/bin"

params.options = [:]
def options    = initOptions(params.options)

process EXTRACT_SPECIES {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"${meta.id}", publish_id:meta.id) }
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-sylph"
        } else {
            conda 'sylph pandas'
        }
    } else {
        conda null
    }
    cache 'lenient'
    scratch true
    
    input:
    tuple val(meta), path(species_obs)
    val(tool)

    output:
    tuple val(meta), stdout, emit: extracted_species
    

    script:
    
    """
    $module_dir/wrangle_file.py $tool $species_obs
    """
}
