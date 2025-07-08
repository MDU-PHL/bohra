// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process STYPE {
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
            conda "${params.conda_path}/bohra-stype"
        } 
        // will need to release stype to conda added in ignore strategy in case people don't use init - at least whole pipeline won't fall down
    } else {
        conda null
    }
    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("typer_${getSoftwareName(task.process)}.txt"), emit: typer
    tuple val(meta), path("version_stype.txt"), emit: version
    

    script:
    """
    echo -e Isolate'\t'Typing_tool'\n'${meta.id}'\t'stype >> tmp.tab
    stype run -c $contigs -px $meta.id
    csvtk csv2tab $meta.id/sistr_filtered.csv > tmp.typer.tab
    paste tmp.tab tmp.typer.tab > raw_typer.tab
    csvtk -t cut -f 'genome,h1,h2,o_antigen,serogroup,serovar,Typing_tool' raw_typer.tab \
    | csvtk -t rename -f 'genome,serogroup,serovar' -n 'Isolate,Serogroup,Serovar' > typer_${getSoftwareName(task.process)}.txt
    echo -e stype'\t'\$CONDA_PREFIX'\t'\$(stype -v) | csvtk add-header -t -n 'tool,conda_env,version' > version_stype.txt
    """
    
}
