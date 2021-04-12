// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

module_dir = moduleDir + "/bin"

params.options = [:]
def options    = initOptions(params.options)
params.species_options = ["Acinetobacter_baumannii", "Campylobacter", "Enterococcus_faecalis", "Enterococcus_faecium", "Escherichia", "Klebsiella", "Salmonella", "Staphylococcus_aureus", "Staphylococcus_pseudintermedius", "Streptococcus_agalactiae", "Streptococcus_pneumoniae", "Streptococcus_pyogenes", "Vibrio_cholerae"]
process AMRFINDER {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    cache 'lenient' 
    cpus options.args2
    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("${meta.id}.out"), emit: amrfinder
    path '*.version.txt'                  , emit: version

    script:
    
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"   
    def organism = params.options.containsValue( meta.species_exp ) ? "--plus --organism $meta.species_exp": ""
    """
    amrfinder -n $contigs -o ${meta.id}.out -t $task.cpus $organism
    echo \$(amrfinder --version 2>&1) | sed -e "s/amrfinder //g" > ${software}.version.txt
    """
    
    
    
}

process COLLATE_AMRFINDER_ISOLATE {
    
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
        
    cache 'lenient'
    input:
    tuple val(meta), path(amrfinder)

    output:
    tuple val(meta), path('summary_matches.txt'), emit: abritamr_matches_isolates
    tuple val(meta), path('summary_partials.txt'), emit: abritamr_partials_isolates
    
    
    script:
    
    """
    ${module_dir}/collate.py $meta.id $amrfinder
    """
    
}


// common process to collate the amr output
process COLLATE_AMRFINDER_MATCHES {
    
    label 'process_medium'
    publishDir ".",
        mode: params.publish_dir_mode
        
    cache 'lenient'
    input:
    val(abritamr_matches_isolates) 
    // val(abritamr_partials_isolates)

    output:
    path 'summary_matches.txt', emit: abritamr_matches
    // path 'summary_partials.csv', emit: abritamr_partials
    // path("${params.prefix}_MMS118.xlsx") optional true
    
    script:
    
    """
    samples=`echo $abritamr_matches_isolates | tr -cd , | wc -c`
    if(( \$samples > 0)); then
    csvtk concat -t \$(echo $abritamr_matches_isolates | sed "s/,/ /g; s/\\[/ /g; s/\\]/ /g") > summary_matches.txt
    else
    cat \$(echo $abritamr_matches_isolates | sed "s/,/\\n/g; s/\\[/ /g; s/\\]/ /g") > summary_matches.txt
    fi
    """
    
}

process COLLATE_AMRFINDER_PARTIALS {
    
    label 'process_medium'
    publishDir ".",
        mode: params.publish_dir_mode
        
    cache 'lenient'
    input:
    // val(abritamr_matches_isolates) 
    val(abritamr_partials_isolates)

    output:
    // path 'summary_matches.csv', emit: abritamr_matches
    path 'summary_partials.txt', emit: abritamr_partials
    // path("${params.prefix}_MMS118.xlsx") optional true
    
    script:
    
    """
    samples=`echo $abritamr_partials_isolates | tr -cd , | wc -c`
    if(( \$samples > 0)); then
    csvtk concat -t \$(echo $abritamr_partials_isolates | sed "s/,/ /g; s/\\[/ /g; s/\\]/ /g") > summary_partials.txt
    else
    cat \$(echo $abritamr_partials_isolates | sed "s/,/\\n/g; s/\\[/ /g; s/\\]/ /g") > summary_partials.txt
    fi
    """
    
}


process MDUIFY_AMRFINDER {
    
    label 'process_medium'
    publishDir ".",
        mode: params.publish_dir_mode
        
    cache 'lenient'
    input:
    tuple path(abritamr_matches), path (abritamr_partials)

    output:
    path("${params.runid}_MMS118.xlsx"), emit: mdu_abritamr
    
    script:
    
    """
    ${module_dir}/collate.py $params.qc $params.db $params.runid $abritamr_matches $abritamr_partials
    """
    
}
