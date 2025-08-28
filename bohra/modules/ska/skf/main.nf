// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SKA_BUILD {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    
    
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/mash" : 'mash') : null) 
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}/${params.conda_prefix}-ska2").exists()) {
            conda "${params.conda_path}/${params.conda_prefix}-ska2"
        } else {
            conda "${moduleDir}/environment.yml"
        }
    } else {
        conda null
    }

    scratch true
    
    cache 'lenient'
    

    input:
    tuple val(meta), val(sequence)

    output:
    tuple val(meta), path("${meta.id}.skf"), emit: skf
    tuple val(meta), path("version_ska.txt"), emit: version

    script:
    input_files = meta.input_type == "pe_reads" ? sequence.join('\t') : sequence
    
    """
    echo -e "$meta.id\t$input_files" > tmp.tsv
    ska build -f tmp.tsv -k $params.ska2_kszise -o $meta.id
    echo -e ska2'\t'\$CONDA_PREFIX'\t'\$(ska --version)'\t'${params.ska_ref} | csvtk add-header -t -n 'tool,conda_env,version,reference' > version_ska.txt
    """
    
    
}
