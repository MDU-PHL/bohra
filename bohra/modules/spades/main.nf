// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SPADES {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    
    cpus options.args2// args2 needs to be cpus for shovill
    cache 'lenient'
    // conda (params.enable_conda ? 'bioconda::shovill=1.1.0' : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container 'https://depot.galaxyproject.org/singularity/fastp:0.20.1--h8b12597_0'
    // } else {
    //     container 'quay.io/biocontainers/fastp:0.20.1--h8b12597_0'
    // }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('contigs.fa'), emit: contigs
    tuple val(meta), path('spades.log'), emit: log

    script:
    """
    tmp_dir=\$(mktemp -d)
    spades.py -1 ${reads[0]} -2 ${reads[1]} -o current -t $task.cpus --isolate --tmp-dir \$tmp_dir
    cp current/scaffolds.fasta contigs.fa
    cp current/spades.log .
    rm -rf \$tmp_dir
    """
    
}
