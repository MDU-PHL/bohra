// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SKESA {
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
    path '*.version.txt'                  , emit: version

    script:
    // Added soft-links to original fastqs for consistent naming in MultiQC
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
    [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
    skesa --fastq ${prefix}_1.fastq.gz,${prefix}_2.fastq.gz --cores $task.cpus --vector_percent 1.0 --contigs_out ${prefix}/contigs.fa
    cp ${prefix}/contigs.fa contigs.fa
    cp ${prefix}/spades.log spades.log
    echo \$(spades.py --version 2>&1) | sed -e "s/spades //g" > ${software}.version.txt
    """
    
}
skesa --fastq "$read1,$read2" \
      --cores "$cpus" \
      --vector_percent 1.0 \
      $opts \
      --contigs_out "$outdir/contigs