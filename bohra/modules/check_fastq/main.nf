// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'


module_dir = moduleDir + "/bin"
params.options = [:]
def options    = initOptions(params.options)

process CHECK_FASTQ {
    tag "$meta.id"
    label 'process_medium'
    
    
    scratch true
    cache 'lenient'
    
    
    input:
    tuple val(meta), path(sequences)

    output:
    tuple val(meta), stdout, emit: pe_check
    
    script:
    
    def read_files = meta.input_type != "pe_reads" ?"$sequences" : "${sequences[0]} ${sequences[1]}"
    """
    $module_dir/check_reads.py $meta.id $read_files
    """
    
}
