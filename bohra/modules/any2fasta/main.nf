// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process ANY2FASTA {
    
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:params.report_outdir, publish_id:params.report_outdir) }


    scratch true
    cache 'lenient'
    
    
    if ( params.enable_conda ) {
        if (file("${params.dependency_prefix}/any2fasta").exists()) {
            conda "${params.dependency_prefix}/any2fasta"
        } 
    } else {
        conda null
    }
    input:
    path(reference)

    output:
    path("reference.fa"), emit: cleaned_reference
    path("version_any2fasta.txt"), emit: version
    

    script:
    """
    any2fasta ${reference} > reference.fa
    echo -e any2fasta'\t'\$CONDA_PREFIX'\t'\$(any2fasta -v 2>&1)'\t'${params.any2fasta_ref} | csvtk add-header -t -n 'tool,conda_env,version,reference' > version_any2fasta.txt
    """
    
}
