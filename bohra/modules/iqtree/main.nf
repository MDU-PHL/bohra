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
        if (file("${params.conda_path}/${params.conda_prefix}-iqtree").exists()) {
            conda "${params.conda_path}/${params.conda_prefix}-iqtree"
        } else {
            conda "${moduleDir}/environment.yml"
    } else {
        conda null
    }

    cache 'lenient'
    scratch true
      
    input:
        path(aln)
        val(full_aln)

    output:
        path('snps.newick'), emit: newick
        path('version_iqtree.txt'), emit: version
    
    script:
    if (full_aln == "no_full_aln") {
        """
        iqtree \\
        -s $aln -pre snps \\
        -m GTR+G4 -bb 1000 -ntmax $task.cpus \\
        -nt AUTO -st DNA
        gotree reroot midpoint -i snps.treefile -o snps.newick
        echo -e iQtree'\t'\$CONDA_PREFIX'\t'\$(iqtree --version | grep version) | csvtk add-header -t -n 'tool,conda_env,version' > version_iqtree.txt
        echo -e gotree'\t'\$CONDA_PREFIX'\t'\$(gotree version)  >> version_iqtree.txt
        """
    } else {
         """
        iqtree -fconst \$(snp-sites -C $full_aln) \\
        -s $aln -pre snps \\
        -m GTR+G4 -bb 1000 -ntmax $task.cpus \\
        -nt AUTO -st DNA
        gotree reroot midpoint -i snps.treefile -o snps.newick
        echo -e iQtree'\t'\$CONDA_PREFIX'\t'\$(iqtree --version | grep version) | csvtk add-header -t -n 'tool,conda_env,version' > version_iqtree.txt
        echo -e gotree'\t'\$CONDA_PREFIX'\t'\$(gotree version)  >> version_iqtree.txt
        """ 
    }
        
    
        
}
