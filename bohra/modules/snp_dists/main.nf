// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SNP_DISTS {
    
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
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
