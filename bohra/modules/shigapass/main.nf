// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'
module_dir = moduleDir + "/bin"
params.options = [:]
def options    = initOptions(params.options)

process SHIGAPASS {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    
    cpus options.args2// args2 needs to be cpus for shovill
    cache 'lenient'
    // errorStrategy 'ignore'
    // scratch true
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/spades" : 'bioconda::spades=3.15.2') : null) 
    if ( params.enable_conda ) {
        if (file("${params.conda_path}/${params.conda_prefix}-shigapass").exists()) {
            conda "${params.conda_path}/${params.conda_prefix}-shigapass"
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
    tuple val(meta), path("version_shigapass.txt"), emit: version
    

    script:
    """
    echo -e Isolate'\t'Typing_tool'\n'${meta.id}'\t'shigapass >> tmp.tab
    echo $contigs > input.tab
    bash ${module_dir}/ShigaPass.sh -l input.tab -o shigapass -t $task.cpus -p ${module_dir}/ShigaPass_DataBases 
    sed 's/;/\t/g' shigapass/ShigaPass_summary.csv > shigapass.tsv
    paste tmp.tab shigapass.tsv > final_output.tsv
    csvtk -t cut -f 'Isolate,Predicted_Serotype,Predicted_FlexSerotype' final_output.tsv > typer_shigapass.txt
    echo -e shigapass'\\t'\$CONDA_PREFIX'\\t'\$(${module_dir}/ShigaPass.sh -v) | csvtk add-header -t -n 'tool,conda_env,version' > version_shigapass.txt
    """
    
}
