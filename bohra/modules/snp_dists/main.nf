// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SNP_DISTS {
    
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:params.report_outdir, publish_id:report_outdir) }
    
    
    scratch true
    if ( params.enable_conda ) {
        if (file("${params.dependency_prefix}/snpdists").exists()) {
            conda "${params.dependency_prefix}/snpdists"
        } 
    } else {
        conda null
    }
    
    cache 'lenient'
    scratch true
    
    input:
    path(core) 

    output:
    path('distances.tsv'), emit: matrix
    path('version_snpdists.txt'), emit: version

    script:
    
    """
    snp-dists $core | csvtk rename -t -f 1 -n Isolate > distances.tsv
    echo -e snp-dists'\t'\$CONDA_PREFIX'\t'\$(snp-dists -v 2>&1)'\t'${params.snpdists_ref} | csvtk add-header -t -n 'tool,conda_env,version,reference' > version_snpdists.txt
    """
    
}
