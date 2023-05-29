// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process QUICKTREE {

    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/quicktree" : 'bioconda::quicktree=2.5') : null) 
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-quicktree"
        } else {
            conda 'bioconda::bioconda::quicktree=2.5 newick_utils'
        }
    } else {
        conda null
    }

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
