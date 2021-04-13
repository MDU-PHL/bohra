// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

module_dir = moduleDir + "/bin"

params.options = [:]
def options    = initOptions(params.options)

process COLLATE_KRAKEN2_ISOLATE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    cache 'lenient' 
    
    input:
    tuple val(meta), path(kraken2)
    
    output:
    tuple val(meta), path('species.txt'), emit: species

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    $module_dir/collate_kraken2.py $meta.id $kraken2 species.txt
    """
    
}

process COLLATE_STATS_ISOLATE {
    
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    cache 'lenient' 
    input:
    tuple val(meta), path(stats), path(seqtk_stats)

    output:
    tuple val(meta), path ("read_assessment.txt"), emit: read_assessment
    
    script:
    """
    ${module_dir}/collate_stats.py $meta.id $stats $seqtk_stats  ${params.reference_fasta} ${params.min_qscore} ${params.min_cov} read_assessment.txt
    """
    
}


process SNIPPY_QC {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    cache 'lenient' 
    
    input:
    tuple val(meta), path(aln)
    
    output:
    tuple val(meta), path('snippy_qc.txt'), emit: snippy_qc

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    $module_dir/snippy_qc.py $meta.id $aln snippy_qc.txt ${params.min_aln}
    """
    
}


process COLLATE_KRAKENS {
    
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode
        
    cache 'lenient'

    input:
    val(krakens) 

    output:
    path "species.txt", emit: collated_species
    
    script:
    def kraken = krakens.join(' ')
    """
    csvtk concat -t -T $kraken > species.txt
    """
}

process COLLATE_SEQDATA {
    
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode
        
    cache 'lenient'
    input:
    val(seqdata) 

    output:
    path "seqdata.txt", emit: collated_seqdata
    
    script:
    def seqs = seqdata.join(' ')
    """
    csvtk concat -t -T $seqs > seqdata.txt
    """

}


process COLLATE_SNIPPY_QCS {
    
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode
        
    cache 'lenient'

    input:
    val(snippy_qcs) 

    output:
    path "core_genome.txt", emit: core_genome
    
    script:
    def sq = snippy_qcs.join(' ')
    """
    csvtk concat -t -T $sq > core_genome.txt
    """
}

process COLLATE_ASM {
    
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    cache 'lenient' 
    input:
    tuple val(meta), path(stats)
    output:
    tuple val(meta), path ("assembly.txt"), emit: assembly
    
    script:
    """
    ${module_dir}/collate_stats.py $meta.id $stats $seqtk_stats  ${params.reference_fasta} ${params.min_qscore} ${params.min_cov} read_assessment.txt
    """


}

process COMPILE {
    
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode
        
    cache 'lenient'
    input:
    val(results) 

    output:
    path "report.html", emit: html
    
    script:
    
    """
    $module_dir/compile.py $params.mode $params.outdir $params.template_dir $launchDir $params.run_kraken -
    """
}