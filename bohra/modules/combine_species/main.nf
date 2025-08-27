// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

module_dir = moduleDir + "/bin"

params.options = [:]
def options    = initOptions(params.options)

process COMBINE_SPECIES_REPORT {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"${meta.id}", publish_id:meta.id) }
    
    cache 'lenient'
    scratch true
    
    input:
    tuple val(meta), val(species_reads), val(species_asm)
    

    output:
    tuple val(meta), path('speciation_report.txt'), emit: species_report
    

    script:
    
    """
    $module_dir/wrangle_file.py $species_reads $species_asm > speciation_report.txt
    """
}
