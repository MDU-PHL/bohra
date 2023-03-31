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
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/seqkit"
        } else {
            conda 'csvtk seqkit=2.1.0 pandas numpy=1.23.1 any2fasta'
        }
    } else {
        conda null
    }
    input:
    tuple val(meta), path(seqkit_stats), path(seqkit_qual), path(reference)

    output:
    tuple val(meta), path ("read_assessment.txt"), emit: read_assessment
    
    script:
    """
    ${module_dir}/collate_stats.py $meta.id $seqkit_stats  \
    $seqkit_qual $launchDir/$params.reference read_assessment.txt
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
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
        
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
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
        
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

process ADD_HEADER_MLST {
    
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    
    cache 'lenient'
    input:
    tuple val(meta), path(mlst)
    
    output:
    tuple val(meta), path('mlst.txt'), emit: mlst

    script:
    
    """
    $module_dir/add_header_mlst.py $mlst 
    """
}


process MOBSUITE_WRANGLE {
    
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    
    cache 'lenient'
    input:
    tuple val(meta), path(mobs)
    
    output:
    tuple val(meta), path('plasmid.txt'), emit: plasmid

    script:
    
    """
    $module_dir/wrangle_mobsuite.py $mobs $meta.id
    """
}


process COLLATE_MOBSUITE {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
        
    cache 'lenient'
    input:
    tuple val(output_name), val(input)

    output:
    path "plasmid.txt", emit: collated_plasmid
    
    script:
    def plasmids = input.join(' ')
    """
    $module_dir/collate_plasmids.py $output_name $plasmids
    """

}

process COLLATE_SNIPPY_QCS {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
        
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
    label 'process_medium'
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    cache 'lenient' 
    input:
    tuple val(meta), path(prokka), path(stats)
    output:
    tuple val(meta), path ("assembly.txt"), emit: assembly
    
    script:
    """
    ${module_dir}/collate_asm.py $meta.id $prokka $stats assembly.txt
    """


}



process COMPILE {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
        
    cache 'lenient'
    input:
    val(results) 

    output:
    path "report.html", emit: html
    
    script:
    def res = results.join(' ')
    """
    $module_dir/compile.py --pipeline $params.mode --launchdir $launchDir \
    --template_dir $params.template_dir  --day ${params.day} --user ${params.user} \
    --isolates $params.isolates --reference $params.reference --job_id $params.job_id \
    --iqtree $params.run_iqtree
    """
}

process COLLATE_ABRITMAR {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
        
    cache 'lenient'
    input:
    tuple val(output_name), val(input)

    output:
    path "${output_name}.txt", emit: collated
    
    script:
    def resistomes = input.join(' ')
    """
    $module_dir/collate_plasmids.py $output_name $resistomes
    """

}