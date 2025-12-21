// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process MASH_TRIANGLE {

    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:params.report_outdir, publish_id:report_outdir) }
    
    
    if ( params.enable_conda ) {
        if (file("${params.dependency_prefix}/mash").exists()) {
            conda "${params.dependency_prefix}/mash"
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
