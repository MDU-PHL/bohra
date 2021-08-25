// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process MASH_TRIANGLE {

    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    
    cache 'lenient'
    scratch true
    
    input:
    val(sketches)

    output:
    path('distances.tab'), emit: mash_distances

    script:
    def input_files = sketches.join(' ')
    """
    mash triangle -C $input_files > distances.tab
    """
        
}
