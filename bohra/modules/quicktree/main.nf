// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process QUICKTREE {

    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:params.report_outdir, publish_id:report_outdir) }
    
    if ( params.enable_conda ) {
        if (file("${params.dependency_prefix}/quicktree").exists()) {
            conda "${params.dependency_prefix}/quicktree"
        } 
    } else {
        conda null
    }

    cache 'lenient'
    scratch true
    input:
    val(distances)

    output:
    path('distance.newick'), emit: newick
    path('version_quicktree.txt'), emit: version

    script:
    """
    quicktree -in m -out t $distances  > tmp.newick
    gotree reroot midpoint -i tmp.newick -o distance.newick
    echo -e quicktree'\t'\$CONDA_PREFIX'\t'\$(quicktree -v)'\t'${params.quicktree_ref} | csvtk add-header -t -n 'tool,conda_env,version,reference' > version_quicktree.txt
    echo -e gotree'\t'\$CONDA_PREFIX'\t'\$(gotree version)'\t'${params.gotree_ref}  >> version_quicktree.txt
    """
        
}
