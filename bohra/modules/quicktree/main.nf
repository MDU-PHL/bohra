// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process QUICKTREE {

    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode
    
    cache 'lenient'
    
    input:
    path(preview_distances)

    output:
    path('preview.newick'), emit: preveiw_tree

    script:
    """
    quicktree -in m -out t $preview_distances  | nw_order -c n - > preview.newick
    """
        
}