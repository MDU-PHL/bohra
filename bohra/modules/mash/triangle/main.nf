// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process MASH_TRIANGLE {

    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/mash" : 'mash') : null) 
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}/${params.conda_prefix}-mash").exists()) {
            conda "${params.conda_path}/${params.conda_prefix}-mash"
        } else {
            conda "${moduleDir}/environment.yml"
        }
    } else {
        conda null
    }

    cache 'lenient'
    scratch true
    
    input:
    val(sketches)

    output:
    path('distances.tsv'), emit: dists

    script:
    def input_files = sketches.join(' ')
    """
    mash triangle -C $input_files > distances.tsv
    """
        
}
