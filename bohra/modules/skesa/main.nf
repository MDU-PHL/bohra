// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SKESA {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"${meta.id}", publish_id:meta.id) }
    
    cache 'lenient'
    // afterScript "rm -fr /tmp/\$USER/*"
    scratch true
    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('contigs.fa'), emit: contigs

    script:
    
    """
    if [ -e $meta.contigs ]
    then
        cp $meta.contigs contigs.fa
    else
        skesa --fastq ${reads[0]},${reads[1]} --cores $task.cpus --vector_percent 1.0 \
        --contigs_out contigs.fa
    fi
    """
    
}
