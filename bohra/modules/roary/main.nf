// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process ROARY {
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"report", publish_id:meta.id) }
    
    cache 'lenient'
    scratch true
    // afterScript "rm -fr /tmp/\$USER/*"
    
    input:
    val(gffs)

    output:
    tuple val(meta), path('summary_statistics.txt'), emit: roary_summary
    tuple val(meta), path('gene_presence_absence.csv'), emit: roary_pres_abs

    script:
    def gffs_str = gffs.join(' ')
    """
    roary -p 36 -f roary $gffs_str
    cp roary/summary_statistics.txt .
    cp roary/gene_presence_absence.csv .
    """
}
