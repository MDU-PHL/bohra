// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process MASH_TRIANGLE {

    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode
    
    cache 'lenient'
    
    input:
    val(sketches)

    output:
    path('preview_distances.tab'), emit: mash_distances

    script:
    // Added soft-links to original fastqs for consistent naming in MultiQC
    def software = getSoftwareName(task.process)
    def input_files = sketches.join(' ')
    """
    mash triangle -C $input_files $params.reference_fasta > preview_distances.tab
    """
        
}
