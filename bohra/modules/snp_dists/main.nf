// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SNP_DISTS {
    
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/snpdists" : 'bioconda::snp-dists=0.8.2 bioconda::csvtk') : null) 
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-snpdists"
        } else {
            conda 'bioconda::snp-dists=0.8.2 bioconda::csvtk'
        }
    } else {
        conda null
    }
    
    cache 'lenient'
    
    input:
    path(core) 

    output:
    path('distances.tab'), emit: distances

    script:
    
    """
    snp-dists $core | csvtk rename -t -f 1 -n Isolate > distances.tab
    """
    
}
