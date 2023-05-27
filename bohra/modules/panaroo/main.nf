// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process PANAROO {
    label 'process_high'
    publishDir "${params.outdir}",
        mode: 'copy',
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"report", publish_id:'report') }
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-panaroo"
        } else {
            conda 'bioconda::panaroo=1.2.9 csvtk'
        }
    } else {
        conda null
    }
    

    cache 'lenient'
    scratch true
    
    
    input:
    val(gffs)

    output:
    path('summary_statistics.txt'), emit: pangenome_summary
    path('gene_presence_absence_roary.csv'), emit: pangenome_csv
    path('panaroo/*')

    script:
    def gffs_str = gffs.join(' ')
    """
    mkdir results
    panaroo -i $gffs_str -o panaroo --clean-mode strict
    csvtk add-header -t -T -n 'Genes,Range,Total' panaroo/summary_statistics.txt > summary_statistics.txt
    cp panaroo/gene_presence_absence_roary.csv .
    ln -sf panaroo $launchDir/
    """
}
