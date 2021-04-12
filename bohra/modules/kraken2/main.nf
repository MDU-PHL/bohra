// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process KRAKEN2 {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    
    cpus options.args2 // args2 needs to be cpus for kraken2
    cache 'lenient'
    // conda (params.enable_conda ? 'bioconda::kraken2=2.1.1' : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container 'https://depot.galaxyproject.org/singularity/fastp:0.20.1--h8b12597_0'
    // } else {
    //     container 'quay.io/biocontainers/fastp:0.20.1--h8b12597_0'
    // }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('kraken2.tab'), emit: kraken2
    path '*.version.txt'                  , emit: version

    script:
    // Added soft-links to original fastqs for consistent naming in MultiQC
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    if (meta.single_end) {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -s $reads ${prefix}.fastq.gz
        kraken2 ${prefix}.trim.fastq.gz \
        --threads $task.cpus \
        --minimum-base-quality 13 \
        --report ${prefix}.kraken2.tab \
        --memory-mapping \
        $options.args
        echo \$(kraken2 --version 2>&1) | sed -e "s/kraken2 //g" > ${software}.version.txt
        """
    } else {
        def fail_fastq = params.save_trimmed_fail ? "--unpaired1 ${prefix}_1.fail.fastq.gz --unpaired2 ${prefix}_2.fail.fastq.gz" : ''
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
        kraken2 --paired ${prefix}_1.fastq.gz  ${prefix}_2.fastq.gz  --threads $task.cpus  --minimum-base-quality 13 --report kraken2.tab --memory-mapping $options.args
        echo \$(kraken2 --version 2>&1) | sed -e "s/kraken2 //g" > ${software}.version.txt
        """
    }
}
