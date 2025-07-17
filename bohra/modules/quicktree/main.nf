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
            conda "${params.conda_path}/bohra-test-quicktree"
        } else {
            conda 'environment.yml'
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
    echo -e quicktree'\t'\$CONDA_PREFIX'\t'\$(quicktree -v) | csvtk add-header -t -n 'tool,conda_env,version' > version_quicktree.txt
    echo -e gotree'\t'\$CONDA_PREFIX'\t'\$(gotree version)  >> version_quicktree.txt
    """
        
}
