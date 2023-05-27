// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process EMMTYPER {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    
    cpus options.args2// args2 needs to be cpus for shovill
    cache 'lenient'
    errorStrategy 'ignore'
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/spades" : 'bioconda::spades=3.15.2') : null) 
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-emmtyper"
        } 
        // will need to release stype to conda added in ignore strategy in case people don't use init - at least whole pipeline won't fall down
    } else {
        conda null
    }
    input:
    tuple val(meta), path(contigs), val(species)

    output:
    tuple val(meta), path('typer.txt'), emit: typer
    

    script:
    """
    echo -e ${meta.id} >> tmp.tab
    emmtpyer $contigs > emmtyper.tab
    paste tmp.tab emmtyper.tab | csvtk add-header -n 'Isolate,Num_clusters,emm_type,emm_like,emm_cluster' > typer.txt
    """
    
}