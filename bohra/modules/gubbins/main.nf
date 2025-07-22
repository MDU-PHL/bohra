// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process GUBBINS {
    
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode
    scratch true
    cache 'lenient'
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/gubbins" : 'gubbins=2.4.1 snp-sites=2.5.1') : null)
    if ( params.enable_conda ) {
        if (file("${params.conda_path}/bohra-gubbins").exists()) {
            conda "${params.conda_path}/bohra-gubbins"
        } else {
            conda "${moduleDir}/environment.yml"
        }
    } else {
        conda null
    }
    input:
    path(cleaned)

    output:
    path('gubbins.aln'), emit: gubbins
    path("version_gubbins.txt"), emit: version
    // tuple val(meta), path('spades.log'), emit: log

    script:
    """
    run_gubbins.py -c $task.cpus  --prefix clean $cleaned
    snp-sites -c clean.filtered_polymorphic_sites.fasta > gubbins.aln
    echo -e gubbins'\t'\$CONDA_PREFIX'\t'\$(run_gubbins.py --version) | csvtk add-header -t -n 'tool,conda_env,version' > version_gubbins.txt
    echo -e snp-sites'\t'\$CONDA_PREFIX'\t'\$(snp-sites -V)  >> version_gubbins.txt
    """
    
}
