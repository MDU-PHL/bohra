// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]

def options    = initOptions(params.options)

process ECTYPER {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    
    cpus options.args2// args2 needs to be cpus for shovill
    cache 'lenient'
    errorStrategy 'ignore'
    scratch true
    
    if ( params.enable_conda ) {
        if (file("${params.dependency_prefix}/ectyper").exists()) {
            conda "${params.dependency_prefix}/ectyper"
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
    tuple val(meta), path("version_ectyper.txt"), emit: version
    

    script:
    """
    echo -e Isolate'\n'${meta.id} >> tmp.tab
    ectyper -i $contigs -o ectyper 
    paste tmp.tab ectyper/output.tsv | csvtk -t cut -f -Name,-Species,-QC > typer_${getSoftwareName(task.process)}.txt
    rm -f tmp.tab
    echo -e ectyper'\t'\$CONDA_PREFIX'\t'\$(ectyper --version)'\t'${params.ectyper_ref} | csvtk add-header -t -n 'tool,conda_env,version,reference' > version_ectyper.txt
    """
    
}
