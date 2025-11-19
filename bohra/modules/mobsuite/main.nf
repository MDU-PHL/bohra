// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process MOBSUITE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    
    cache 'lenient'
    scratch true
    
    if ( params.enable_conda ) {
        if (file("${params.dependency_prefix}/mob_suite").exists()) {
            conda "${params.dependency_prefix}/mob_suite"
        } else {
            conda "${moduleDir}/environment.yml"
        }
    } else {
        conda null
    }

    input:
    tuple val(meta), path(asm)

    output:
    tuple val(meta), path('contig_report.txt'), emit: contig_report
    tuple val(meta), path('mobtyper_results.txt'), emit: mobs
    tuple val(meta), path('*.fasta'), optional: true
    tuple val(meta), path("version_mobsuite.txt"), emit: version
    // tuple val(meta), path('spades.log'), emit: log

    script:
    def db = params.mobsuite_db != 'no_db'  ?  "--database $params.mobsuite_db" : ""
    """
    mob_recon -i $asm -s $meta.id -n $task.cpus -o mob $db
    cp mob/* .
    if [ ! -f mobtyper_results.txt ];then
        touch mobtyper_results.txt
    fi
    echo -e mobsuite'\t'\$CONDA_PREFIX'\t'\$(mob_recon -V)'\t'$db'\t'${params.mobsuite_ref} | csvtk add-header -t -n 'tool,conda_env,version,database,reference' > version_mobsuite.txt
    """
    
}
