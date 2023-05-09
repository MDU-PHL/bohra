// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SNIPPY_CORE {
    
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode
    
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/snippy" : 'bioconda::snippy=4.4.5') : null) 
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-snippy"
        } else {
            conda 'bioconda::snippy=4.4.5'
        }
    } else {
        conda null
    }
    cpus options.args2// args2 needs to be cpus for shovill
    cache 'lenient'
    
    input:
    val(alns) // this needs to be a list of sample! not .aln since snippy core uses relative path and the name of the folder to name results!
    path(reference)
    output:
    path('core.aln'), emit: core_aln
    path('core.full.aln'), emit: core_full_aln
    path('core.vcf'), emit: core_vcf
    path('core.txt'), emit: core_stats

    script:
    
    def mask_string = params.mask_string != 'no_mask' ? "--mask $launchDir/${params.mask_string}" : ""
    def core = alns.join(' ')
    """
    snippy-core --ref $reference ${mask_string} $core
    """
    
}
