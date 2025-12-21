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
    
    if ( params.enable_conda ) {
        if (file("${params.dependency_prefix}/gubbins").exists()) {
            conda "${params.dependency_prefix}/gubbins"
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
    echo -e gubbins'\t'\$CONDA_PREFIX'\t'\$(run_gubbins.py --version)'\t'${params.gubbins_ref} | csvtk add-header -t -n 'tool,conda_env,version,reference' > version_gubbins.txt
    echo -e snp-sites'\t'\$CONDA_PREFIX'\t'\$(snp-sites -V)  >> version_gubbins.txt
    """
    
}
