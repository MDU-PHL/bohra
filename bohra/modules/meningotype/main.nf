// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process MENINGOTYPE {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    
    
    cache 'lenient'
    scratch true
    errorStrategy 'ignore'
    
    if ( params.enable_conda ) {
        if (file("${params.dependency_prefix}/meningotype").exists()) {
            conda "${params.dependency_prefix}/meningotype"
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
    tuple val(meta), path("version_meningotype.txt"), emit: version

    script:
    """
    echo -e Isolate'\n'${meta.id} >> tmp.tab
    meningotype --all $contigs  > meningotype.tab
    paste tmp.tab meningotype.tab | csvtk -t rename -f SEROGROUP -n Serogroup | csvtk -t cut -f -SAMPLE_ID,-MLST >typer_${getSoftwareName(task.process)}.txt
    rm -f tmp.tab
    echo -e meningotype'\t'\$CONDA_PREFIX'\t'\$(meningotype --version)'\t'${params.meningotype_ref} | csvtk add-header -t -n 'tool,conda_env,version,reference' > version_meningotype.txt
    """
    
}
