// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

module_dir = moduleDir + "/bin"

params.options = [:]
def options    = initOptions(params.options)


process ABRITAMR {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"${meta.id}/abritamr/current", publish_id:meta.id) }
    // publishDir "abritamr",
    //     mode: 'link',
    //     saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }

    cache 'lenient' 
    scratch true
    // afterScript "rm -fr /tmp/\$USER/*"
    
    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("summary_virulence.txt"), emit: abritamr_virulence
    tuple val(meta),path('summary_matches.txt'), emit: abritamr_matches
    tuple val(meta),path('summary_partials.txt'), emit: abritamr_partials

    script:
    
    def organism = params.species_options.any { it.contains(meta.species_exp) } ? "-sp $meta.species_exp": ""
    """
    abritamr run -c $contigs -px ${meta.id} -j $task.cpus $organism
    cp ${meta.id}/*.txt .
    """
}

process COMBINE_AMR {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    
    cache 'lenient' 
    scratch true

    input:
    tuple val(meta), path(matches), path(partials)

    output:
    tuple val(meta), path("resistome.txt"), emit: resistome
    
    script:
    """
    $module_dir/combine.py ${meta.id} $matches $partials
    """

}

process COLLATE_ABRITAMR {
    
    label 'process_low'
        
    cache 'lenient'
    input:
    val(abritamr_matches_isolates) 
    
    output:
    path 'summary_matches.txt', emit: abritamr_matches
    
    
    script:
    
    """
    $module_dir/concat.py summary_matches.txt $abritamr_matches_isolates
    """
    
}

process COLLATE_ABRITAMR_PARTIALS {
    
    label 'process_low'
      
    cache 'lenient'
    input:
    val(abritamr_partials_isolates)

    output:
    path 'summary_partials.txt', emit: abritamr_partials
    
    script:
    
    """
    $module_dir/concat.py summary_partials.txt $abritamr_partials_isolates
    """
    
}


