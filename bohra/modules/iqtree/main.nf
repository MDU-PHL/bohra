// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)


process IQTREE {

    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/iqtree2" : 'iqtree=2.1.4 snp-sites=2.5.1') : null)
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-iqtree"
        } else {
            conda 'iqtree=2.1.4 snp-sites=2.5.1'
        }
    } else {
        conda null
    }

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
