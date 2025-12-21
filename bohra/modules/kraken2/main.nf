// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

module_dir = moduleDir + "/bin"

params.options = [:]
def options    = initOptions(params.options)

process KRAKEN2 {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"${meta.id}", publish_id:meta.id) }
    
    if ( params.enable_conda ) {
        if (file("${params.dependency_prefix}/kraken2").exists()) {
            conda "${params.dependency_prefix}/kraken2"
        } 
    } else {
        conda null
    }
    cache 'lenient'
    scratch true
    
    input:
    tuple val(meta), path(sequences)
    

    output:
    tuple val(meta), path('kraken2.tab'), emit: species_raw
    tuple val(meta), path('species.txt'), emit: species
    tuple val(meta), path('version_kraken2.txt'), emit: version

    script:
    def input_file = meta.input_type != "pe_reads" ?"$sequences" : "--paired ${sequences[0]} ${sequences[1]}"
    """
    kraken2  $input_file \
    --threads $task.cpus \
    --report kraken2.tab \
    --output - \
    -db $params.kraken2_db \
    --memory-mapping 2> /dev/null

    $module_dir/collate_kraken2.py $meta.id kraken2.tab species.txt
    echo -e kraken2'\t'\$CONDA_PREFIX'\t'\$(kraken2 --version | grep version)'\t'$params.kraken2_db'\t'${params.kraken2_ref} | csvtk add-header -t -n 'tool,conda_env,version,database,reference' > version_kraken2.txt
    """
}
