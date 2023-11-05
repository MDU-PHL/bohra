// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process NGMASTER {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    
    cpus options.args2// args2 needs to be cpus for shovill
    cache 'lenient'
    errorStrategy 'ignore'
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-ngmaster"
        } else {
            conda 'ngmaster csvtk'
        }
       
    } else {
        conda null
    }
    input:
    tuple val(meta), path(contigs), val(species)

    output:
    tuple val(meta), path("typer_${getSoftwareName(task.process)}.txt"), emit: typer
    
    

    script:
    """
    echo -e Isolate'\n'${meta.id} >> tmp.tab
    ngmaster  $contigs  > ngmaster.tab
    paste tmp.tab ngmaster.tab | csvtk -t cut -f -2> typer_${getSoftwareName(task.process)}.txt
    rm -f tmp.tab
    """
    
}
