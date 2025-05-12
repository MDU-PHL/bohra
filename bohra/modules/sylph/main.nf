// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

module_dir = moduleDir + "/bin"

params.options = [:]
def options    = initOptions(params.options)

process SYLPH {
    tag "$meta.id"
    label 'process_medium'
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
    // scratch true
    
    input:
    tuple val(meta), path(sequences)

    output:
    tuple val(meta), path('sylph_raw.tsv'), emit: species_raw
    tuple val(meta), path('sylph.tsv'), emit: species



    script:
    def input = meta.input_type == "pe_reads" ?  "-1 ${sequences[0]} -2 ${sequences[1]}":  "$sequences"
    """
    sylph profile $params.sylph_db $input -o sylph_raw.tsv  -t $task.cpus 
    $module_dir/wrangle_file.py $meta.id sylph_raw.tsv > sylph.tsv
    """
}
