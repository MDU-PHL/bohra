// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SKA_BUILD {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    
    
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/mash" : 'mash') : null) 
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-ska2"
        } else {
            conda 'ska2'
        }
    } else {
        conda null
    }

    scratch true
    
    cache 'lenient'
    

    input:
    tuple val(meta), path(sequence)

    output:
    tuple val(meta), path("${meta.id}.skf"), emit: skf

    script:
    input_files = meta.input_type == "pe_reads" ? "$sequence[0]\t$sequence[1]" : sequence
    """
    echo -e "$meta.id\t$input_files" > tmp.tsv
    ska build -f tmp.csv -k $params.ska2_kszise -o $meta.id
    """
    
}
