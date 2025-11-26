// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process INSERTIQR {
    tag "$meta.id"
    label 'process_high'
     publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"${meta.id}", publish_id:meta.id) }  
    cache 'lenient'
    
    if ( params.enable_conda ) {
        if (file("${params.dependency_prefix}/snippy").exists()) {
            conda "${params.dependency_prefix}/snippy"
        } 
    } else {
        conda null
    }
    scratch true
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(asm), val(reads)

    output:    
    tuple val(meta), path('insertiqr.txt'),  emit: stats

    script:
    """
    snippy --outdir ${meta.id} --ref $asm --R1 ${reads[0]} --R2 ${reads[1]} --cpus $task.cpus --force 
    samtools depth -aa ${meta.id}/snps.bam > depth.txt
    x=\$(csvtk summary -t -H -f 3:q1,3:q3 depth.txt | csvtk -t -H mutate2 -e '\$2-\$1' | cut -f3)
    y=\$(samtools stats ${meta.id}/snps.bam | grep 'insert size average' | cut -f 3)
    echo ${meta.id},\$x,\$y >> tmp.csv
    csvtk add-header -n 'Isolate,IQR depth,Insert size' tmp.csv | csvtk csv2tab > insertiqr.txt
    """
    
}

