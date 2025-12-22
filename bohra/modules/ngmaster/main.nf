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
    
    cache 'lenient'
    errorStrategy 'ignore'
    scratch true
    if ( params.enable_conda ) {
        if (file("${params.dependency_prefix}/ngmaster").exists()) {
            conda "${params.dependency_prefix}/ngmaster"
        } 
       
    } else {
        conda null
    }
    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("typer_${getSoftwareName(task.process)}.txt"), emit: typer
    tuple val(meta), path("version_ngmaster.txt"), emit: version
    

    script:
    """
    echo -e Isolate'\n'${meta.id} >> tmp.tab
    ngmaster  $contigs  > ngmaster.tab
    paste tmp.tab ngmaster.tab | csvtk -t cut -f -2> typer_${getSoftwareName(task.process)}.txt
    rm -f tmp.tab
    echo -e ngmaster'\t'\$CONDA_PREFIX'\t'\$(ngmaster --version)'\t'${params.ngmaster_ref} | csvtk add-header -t -n 'tool,conda_env,version,reference' > version_ngmaster.txt
    """
    
}
