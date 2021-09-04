// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process PANAROO {
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"report", publish_id:'report') }
    
    cache 'lenient'
    scratch true
    // afterScript "rm -fr /tmp/\$USER/*"
    
    input:
    val(gffs)

    output:
    path('summary_statistics.txt'), emit: roary_summary
    path('gene_presence_absence.csv'), emit: roary_csv

    script:
    def gffs_str = gffs.join(' ')
    """
    mkdir results
    panaroo -i *.gff -o results --clean-mode strict
    csvtk add-header -t -T -n 'Genes,Range,Total' roary/summary_statistics.txt > summary_statistics.txt
    cp roary/gene_presence_absence.csv .
    """
}
