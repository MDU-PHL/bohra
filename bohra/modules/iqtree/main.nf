// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

def module_dir = moduleDir + "/bin"


process IQTREE {

    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode
    
    cache 'lenient'
        
    input:
        path(aln)

    output:
        path('core.newick'), emit: newick

    script:    
    """
    $module_dir/iqtree_generator.sh ${params.reference_fasta} $aln core $task.cpus
    cp core.treefile core.newick
    """
        
}
