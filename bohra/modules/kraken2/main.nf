// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process KRAKEN2 {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"${meta.id}", publish_id:meta.id) }
    
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/kraken2" : 'kraken2=2.1.2') : null) 
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-kraken2"
        } else {
            conda 'kraken2=2.1.2'
        }
    } else {
        conda null
    }
    cache 'lenient'
    scratch true
    // afterScript "rm -fr /tmp/\$USER/*"
    
    input:
    tuple val(meta), path(reads), path(kraken)

    output:
    tuple val(meta), path('kraken2.tab'), emit: kraken2


    script:
    def read_files = meta.single_end ?"$reads" : "--paired ${reads[0]} ${reads[1]}"
    """
    kraken2  $read_files \
    --threads $task.cpus \
    --report kraken2.tab \
    --output - \
    -db $kraken \
    --memory-mapping \
    $options.args 2> /dev/null
    """
}
