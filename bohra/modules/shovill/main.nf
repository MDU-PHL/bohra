// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SHOVILL {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"${meta.id}", publish_id:meta.id) }
    
    
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/shovill" : 'shovill=1.1.0') : null) 
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-shovill"
        } else {
            conda 'bioconda::shovill=1.1.0'
        }
    } else {
        conda null
    }

    cache 'lenient'
    // afterScript "rm -fr /tmp/\$USER/*"
    scratch true
    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('contigs.fa') , emit: contigs

    script:
    
    """
    if [ -e $meta.contigs ]
    then
        cp $meta.contigs contigs.fa
    else
        shovill --R1 ${reads[0]} --R2 ${reads[1]} --outdir current --cpus $task.cpus --ram 16 
        cp current/contigs.fa contigs.fa
    fi
    """
    
}
