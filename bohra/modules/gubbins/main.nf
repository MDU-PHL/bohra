// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process GUBBINS {
    
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode
    
    cache 'lenient'
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/gubbins" : 'gubbins=2.4.1 snp-sites=2.5.1') : null)
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-gubbins"
        } else {
            conda 'gubbins=2.4.1 snp-sites=2.5.1 csvtk'
        }
    } else {
        conda null
    }
    input:
    path(cleaned)

    output:
    path('gubbins.aln'), emit: gubbins

    // tuple val(meta), path('spades.log'), emit: log

    script:
    """
    run_gubbins.py -c $task.cpus  --prefix clean $cleaned
    snp-sites -c clean.filtered_polymorphic_sites.fasta > gubbins.aln
    """
    
}
