// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'
module_dir = moduleDir + "/bin"
params.options = [:]
def options    = initOptions(params.options)

process SONNEITYPE {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    
    cpus options.args2// args2 needs to be cpus for shovill
    cache 'lenient'
    errorStrategy 'ignore'
    // scratch true
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/spades" : 'bioconda::spades=3.15.2') : null) 
    if ( params.enable_conda ) {
        if (file("${params.conda_path}/${params.conda_prefix}-sonneitype").exists()) {
            conda "${params.conda_path}/${params.conda_prefix}-sonneitype"
        } else {
            conda "${moduleDir}/environment.yml"
        }
        
    } else {
        conda null
    }
    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("typer_${getSoftwareName(task.process)}.txt"), emit: typer
    tuple val(meta), path("sonneitype.txt"), emit: raw
    tuple val(meta), path("version_sonneitype_mykrobe.txt"), emit: version
    

    script:
    """
    mkdir ${meta.id}
    mykrobe predict --sample ${meta.id} --species sonnei --format json --out ${meta.id}/${meta.id}.json --seq ${reads[0]} ${reads[1]}
    ${module_dir}/parse_mykrobe_predict.py --jsons ${meta.id}/${meta.id}.json --alleles ${module_dir}/alleles.txt --prefix ${meta.id}/sonneitype
    csvtk -t rename -f 'genome,final genotype,name' ${meta.id}/sonneitype_predictResults.tsv > typer_sonneitype.txt
    cp ${meta.id}/sonneitype_predictResults.tsv sonneitype.txt
    echo -e mykrobe (sonneitype)'\t'\$CONDA_PREFIX'\t'\$(mykrobe --version) | csvtk add-header -t -n 'tool,conda_env,version' > version_mykrobe.txt
    """
    
}
