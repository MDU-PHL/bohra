// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process AURICLASS {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    
    cpus options.args2// args2 needs to be cpus for shovill
    cache 'lenient'
    scratch true
    
    if ( params.enable_conda ) {
        if (file("${params.dependency_prefix}/auriclass").exists()) {
            conda "${params.dependency_prefix}/auriclass"
        } 
        
    } else {
        conda null
    }
    input:
    tuple val(meta), path(sequences)

    output:
    tuple val(meta), path("typer_${getSoftwareName(task.process)}.txt"), emit: typer
    tuple val(meta), path("version_auriclass.txt"), emit: version


    script:
    def reads = sequences.join(' ')
    """
    auriclass --name  ${meta.id} -o tmp.tsv $reads
    csvtk -t rename -f Sample -n Isolate tmp.tsv > typer_${getSoftwareName(task.process)}.txt
    echo -e auriclass'\t'\$CONDA_PREFIX'\t'\$(auriclass --version)'\t'${params.sccme_ref} | csvtk add-header -t -n 'tool,conda_env,version,reference' > version_auriclass.txt
    """
    
}
