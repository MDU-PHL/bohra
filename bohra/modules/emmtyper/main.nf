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
    scratch true
    
    if ( params.enable_conda ) {
        if (file("${params.dependency_prefix}/emmtyper").exists()) {
            conda "${params.dependency_prefix}/emmtyper"
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
    tuple val(meta), path("version_emmtyper.txt"), emit: version
    
    script:
    """
    echo -e ${meta.id} >> tmp.tab
    emmtyper $contigs > emmtyper.tab
    paste tmp.tab emmtyper.tab | csvtk -t cut -f -2 |csvtk -t add-header -n 'Isolate,Num_clusters,emm_type,emm_like,emm_cluster' > typer_${getSoftwareName(task.process)}.txt
    rm -f tmp.tab
    echo -e emmtyper'\t'\$CONDA_PREFIX'\t'\$(emmtyper --version)'\t'${params.emmtyper_ref} | csvtk add-header -t -n 'tool,conda_env,version,reference' > version_emmtyper.txt
    """
    
}
