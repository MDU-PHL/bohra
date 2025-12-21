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
        if (file("${params.dependency_prefix}/coresnpfilter").exists()) {
            conda "${params.dependency_prefix}/coresnpfilter"
        } 
    } else {
        conda null
    }
    input:
    path(core_full_aln)

    output:
    // tuple val(meta), path("${meta.id}/*"), emit: snippy_dir
    path("core.filtered.aln"), emit: aln
    path("core_snp_table.tsv"), emit: table
    path("version_coresnp.txt"), emit: version
    script:
    """
    coresnpfilter -c ${params.fuzzy_core_prop} --table core_snp_table.tsv -e $core_full_aln > core.filtered.aln 
    echo -e core-snp-filter'\t'\$CONDA_PREFIX'\t'\$(coresnpfilter --version)'\t'${params.core_snp_filter_ref} | csvtk add-header -t -n 'tool,conda_env,version,reference' > version_coresnp.txt
    """
    
}
