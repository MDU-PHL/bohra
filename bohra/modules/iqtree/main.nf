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
        path(full_aln)

    output:
        path('core.newick'), emit: newick

    script:    
    """
    iqtree -fconst \$(snp-sites -C $full_aln) \\
    -s $aln -pre core \\
    -m GTR+G4 -bb 1000 -ntmax $task.cpus \\
    -nt AUTO -st DNA
    cp core.treefile core.newick
    """
        
}
