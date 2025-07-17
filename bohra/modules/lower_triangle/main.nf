// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

module_dir = moduleDir + "/bin"

params.options = [:]
def options    = initOptions(params.options)

process LOWER_TRIANGLE {
    label 'process_high'
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
