// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)


process VERYFASTTREE {

    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/iqtree2" : 'iqtree=2.1.4 snp-sites=2.5.1') : null)
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-veryfasttree"
        } else {
            conda 'veryfasttree'
        }
    } else {
        conda null
    }

    cache 'lenient'
    scratch true
      
    input:
        path(aln)
        // path(full_aln)

    output:
        path('snps.newick'), emit: newick
        path('version_veryfasttree.txt'), emit: version
        
    script:    
    """
    VeryFastTree -nt -gamma -gtr -threads $task.cpus $aln > tmp.newick
    gotree reroot midpoint -i tmp.newick -o snps.newick
    echo -e VeryFastTree'\t'\$CONDA_PREFIX'\t'\$(VeryFastTree --help | head -n 1 2>&1) | csvtk add-header -t -n 'tool,conda_env,version' > version_veryfasttree.txt
    echo -e gotree'\t'\$CONDA_PREFIX'\t'\$(gotree version)  >> version_veryfasttree.txt
    """
        
}
