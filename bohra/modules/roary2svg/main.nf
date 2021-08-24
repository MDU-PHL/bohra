// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

module_dir = moduleDir + "/bin"

process ROARY2SVG {
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"report", publish_id:'report') }
    
    cache 'lenient'
    scratch true
    // afterScript "rm -fr /tmp/\$USER/*"
    
    input:
    path(csv)

    output:
    path('pan_genome.svg'), emit: pan_genome

    script:
    
    """
    $module_dir/roary2svg.pl $csv > pan_genome.svg
    """
}