// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

def module_dir = moduleDir + "/bin"


process IQTREE {

    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    cache 'lenient'
    scratch true
      
    input:
        path(aln)
        path(reference)

    output:
        path('core.newick'), emit: newick

    script:    
    """
    $module_dir/iqtree_generator.sh $reference $aln core $task.cpus
    cp core.treefile core.newick
    """
        
}
