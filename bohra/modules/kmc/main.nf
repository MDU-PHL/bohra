// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process KMC {
    tag "$meta.id"
    label 'process_high'
    // afterScript "rm -fr /tmp/\$USER/*"
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-kmc"
        } else {
            conda 'kmc csvtk'
        }
    } else {
        conda null
    }

    scratch true
    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('est_genome_size_kmer.txt'), emit: genome_size


    script:
    
    """
    tmp_dir=\$(mktemp -d)
    kmc  -t$task.cpus $options.args ${reads[0]} \$tmp_dir/kmc \$tmp_dir | grep 'No. of unique counted k-mers' | cut -f2 -d: > est_genome_size_kmer.txt
    rm -rf \$tmp_dir
    """
    
}
