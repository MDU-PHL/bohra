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
    scratch true
    if ( params.enable_conda ) {
        if (file("${params.dependency_prefix}/lissero").exists()) {
            conda "${params.dependency_prefix}/lissero"
        } else {
            conda "${moduleDir}/environment.yml"
        }
        
    } else {
        conda null
    }
    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("typer_${getSoftwareName(task.process)}.txt"), emit: typer
    tuple val(meta), path("version_lissero.txt"), emit: version
    script:
    """
    echo -e Isolate'\t'Typer_toolname'\n'${meta.id}'\t'${getSoftwareName(task.process)} >> tmp.tab
    lissero $contigs | sed 's/contigs\\.fa/$meta.id/g'  > lissero.tab
    paste tmp.tab lissero.tab | csvtk -t rename -f SEROTYPE -n Serotype | csvtk -t cut -f -ID,-COMMENT > typer_${getSoftwareName(task.process)}.txt
    rm -f tmp.tab
    echo -e lissero'\t'\$CONDA_PREFIX'\t'\$(lissero --version)'\t'${params.lissero_ref} | csvtk add-header -t -n 'tool,conda_env,version,reference' > version_lissero.txt
    """
    
}
