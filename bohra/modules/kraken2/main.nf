// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process KRAKEN2 {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    
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
    
    script:
    if (meta.single_end) {
        """
        kraken2 $reads \
        --threads $task.cpus \
        --minimum-base-quality 13 \
        --report ${prefix}.kraken2.tab \
        --memory-mapping \
        $options.args
        """
    } else {
        def fail_fastq = params.save_trimmed_fail ? "--unpaired1 ${reads[0]} --unpaired2 ${reads[1]}" : ''
        """
        kraken2 --paired ${reads[0]}  ${reads[1]}  --threads $task.cpus  --minimum-base-quality 13 --report kraken2.tab --memory-mapping $options.args
        """
    }
}
