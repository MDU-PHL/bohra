// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process CORE_SNP_FILTER {
    
    label 'process_medium'
    publishDir "${params.outdir}/report",
        mode: params.publish_dir_mode
        
    
    // scratch true
    cache 'lenient'
    
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/snippy" : 'bioconda::snippy=4.4.5') : null) 
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-core-snp-filter"
        } else {
            conda 'bioconda::core-snp-filter'
        }
    } else {
        conda null
    }
    input:
    path(core_full_aln)

    output:
    // tuple val(meta), path("${meta.id}/*"), emit: snippy_dir
    path("core.filtered.aln"), emit: aln
        
    script:
    """
    coresnpfilter -c ${params.coresnpfilter_prop} -e $core_full_aln > core.filtered.aln 
    echo -e core-snp-filter'\t'\$CONDA_PREFIX'\t'\$(coresnpfilter --version) | csvtk add-header -t -n 'tool,conda_env,version' > version_spades.txt
    """
    
}
