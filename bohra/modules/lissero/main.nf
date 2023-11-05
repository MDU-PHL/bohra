// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process LISSERO {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    
    cpus options.args2// args2 needs to be cpus for shovill
    cache 'lenient'
    // errorStrategy 'ignore'
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/spades" : 'bioconda::spades=3.15.2') : null) 
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-lissero"
        } else {
            conda 'lissero csvtk'
        }
        
    } else {
        conda null
    }
    input:
    tuple val(meta), path(contigs), val(species_obs)

    output:
    tuple val(meta), path("typer_${getSoftwareName(task.process)}.txt"), emit: typer

    script:
    """
    echo -e Isolate'\n'${meta.id} >> tmp.tab
    lissero $contigs | sed 's/contigs\\.fa/$meta.id/g'  > lissero.tab
    paste tmp.tab lissero.tab | csvtk -t rename -f SEROTYPE -n Serotype | csvtk -t cut -f -ID,-COMMENT > typer_${getSoftwareName(task.process)}.txt
    rm -f tmp.tab
    """
    
}
