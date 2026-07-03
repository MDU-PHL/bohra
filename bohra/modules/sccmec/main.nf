// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SCCMEC {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    
    cpus options.args2// args2 needs to be cpus for shovill
    cache 'lenient'
    scratch true
    
    if ( params.enable_conda ) {
        if (file("${params.dependency_prefix}/sccmec").exists()) {
            conda "${params.dependency_prefix}/sccmec"
        } 
        
    } else {
        conda null
    }
    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("typer_${getSoftwareName(task.process)}.txt"), emit: typer
    tuple val(meta), path("${meta.id}.tsv"), emit: rawtyper
    tuple val(meta), path("version_sccme.txt"), emit: version
    script:
    """
    sccmec --input $contigs --prefix ${meta.id}
    csvtk -t rename -f sample -n Isolate ${meta.id}.tsv | csvtk -t cut -f Isolate,type,subtype,mecA,comment > typer_${getSoftwareName(task.process)}.txt
    echo -e sccme'\t'\$CONDA_PREFIX'\t'\$(sccme -v 2>&1)'\t'${params.sccme_ref} | csvtk add-header -t -n 'tool,conda_env,version,reference' > version_sccme.txt
    """
    
}
