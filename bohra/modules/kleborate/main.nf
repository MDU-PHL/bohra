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
    
    cache 'lenient'
    scratch true
    errorStrategy 'ignore'
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/spades" : 'bioconda::spades=3.15.2') : null) 
    if ( params.enable_conda ) {
        if (file("${params.dependency_prefix}/kleborate").exists()) {
            conda "${params.dependency_prefix}/kleborate"
        }
        
    } else {
        conda null
    }
    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("kleborate_raw.tab"), emit: raw_klebs
    tuple val(meta), path("typer_${getSoftwareName(task.process)}.txt"), emit: typer
    tuple val(meta), path("version_kleborate.txt"), emit: version
    

    script:
    """
    echo -e Isolate'\t'Typing_tool'\n'${meta.id}'\t'kleborate >> tmp.tab
    kleborate -k -o kleborate_raw.tab -a $contigs
    paste tmp.tab kleborate_raw.tab | csvtk -t rename -f species -n Species | csvtk -t cut -f Isolate,Typing_tool,RmpADC,RmST,rmpA2,wzi,K_locus,K_type,O_locus,O_type > typer_${getSoftwareName(task.process)}.txt
    rm -f tmp.tab
    echo -e kleborate'\t'\$CONDA_PREFIX'\t'\$(kleborate --version)'\t'${params.kleborate_ref} | csvtk add-header -t -n 'tool,conda_env,version,reference' > version_kleborate.txt
    """
    
}

