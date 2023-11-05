// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process KLEBORATE {
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
            conda "${params.conda_path}/bohra-kleborate"
        } else {
            conda 'bioconda::kleborate bioconda::csvtk'
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
    kleborate -o kleborate.tab -a $contigs
    paste tmp.tab kleborate.tab | csvtk -t rename -f species -n Species | csvtk -t cut -f -strain,-ST > typer_${getSoftwareName(task.process)}.txt
    rm -f tmp.tab
    """
    
}
