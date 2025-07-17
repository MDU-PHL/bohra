// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process CORE_SNP_FILTER {
    
    label 'process_medium'
    publishDir "${params.outdir}/report",
        mode: params.publish_dir_mode
        
    
    scratch true
    cache 'lenient'
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}/bohra-core-snp-filter").exists()) {
            conda "${params.conda_path}/bohra-core-snp-filter"
        } else {
            conda 'environment.yml'
        }
    } else {
        conda null
    }
    input:
    path(core_full_aln)

    output:
    // tuple val(meta), path("${meta.id}/*"), emit: snippy_dir
    path("core.filtered.aln"), emit: aln
    path("version_coresnp.txt"), emit: version
    script:
    """
    coresnpfilter -c ${params.fuzzy_core_prop} -e $core_full_aln > core.filtered.aln 
    echo -e core-snp-filter'\t'\$CONDA_PREFIX'\t'\$(coresnpfilter --version) | csvtk add-header -t -n 'tool,conda_env,version' > version_coresnp.txt
    """
    
}
