// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'
module_dir = moduleDir + "/bin"
params.options = [:]
def options    = initOptions(params.options)



process EXTRACT_SPECIES {
    tag "$meta.id"
    label 'process_low'

    input:
        tuple val(meta), path(kraken2)
    output:
        tuple val(meta), env(species), emit: species_obs

    script:
    """
    species=\$($module_dir/extract_species.py  $kraken2)
    """

}
