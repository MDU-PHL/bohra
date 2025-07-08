// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

module_dir = moduleDir + "/bin"

params.options = [:]
def options    = initOptions(params.options)

process LOWER_TRIANGLE {
    label 'process_high'
    // publishDir "${params.outdir}",
    //     mode: params.publish_dir_mode,
    //     saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"report", publish_id:meta.id) }
    
    // if ( params.enable_conda ) {
    //     if (file("${params.conda_path}").exists()) {
    //         conda "${params.conda_path}/bohra-sylph"
    //     } else {
    //         conda 'sylph pandas'
    //     }
    // } else {
    //     conda null
    // }
    cache 'lenient'
    scratch true
    
    input:
    path(distances)

    

    output:
    path('distance_triangle.txt'), emit: triangle
    

    script:
    
    """
    $module_dir/lower_triangle.py $distances > distance_triangle.txt
    """
}
