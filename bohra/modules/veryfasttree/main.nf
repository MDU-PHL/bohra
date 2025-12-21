// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)


process VERYFASTTREE {

    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:params.report_outdir, publish_id:report_outdir) }
    
    
    if ( params.enable_conda ) {
        if (file("${params.dependency_prefix}/veryfasttree").exists()) {
            conda "${params.dependency_prefix}/veryfasttree"
        } else {
            conda "${moduleDir}/environment.yml"
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
    echo -e VeryFastTree'\t'\$CONDA_PREFIX'\t'\$(VeryFastTree --help | head -n 1 2>&1)'\t'${params.veryfasttree_ref} | csvtk add-header -t -n 'tool,conda_env,version,reference' > version_veryfasttree.txt
    echo -e gotree'\t'\$CONDA_PREFIX'\t'\$(gotree version)'\t'${params.gotree_ref}  >> version_veryfasttree.txt
    """
        
}
