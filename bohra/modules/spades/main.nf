// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SPADES {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }

    cache 'lenient'
    
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/spades" : 'bioconda::spades=3.15.2') : null) 
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-spades"
        } else {
            conda 'bioconda::spades=3.13 csvtk python=3.7'
        }
    } else {
        conda null
    }
    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('contigs.fa'), emit: contigs
    // tuple val(meta), path('spades.log'), emit: log

    script:
    """
    if [ -e $meta.contigs ]
    then
        cp $meta.contigs contigs.fa
    else
        tmp_dir=\$(mktemp -d)
        spades.py -1 ${reads[0]} -2 ${reads[1]} -o current -t $task.cpus $options.args2 --tmp-dir \$tmp_dir
        cp current/contigs.fasta contigs.fa
        rm -rf \$tmp_dir
    fi
    """
    
}
