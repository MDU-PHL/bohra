// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process PANAROO {
    label 'process_high'
    publishDir "${params.outdir}",
        mode: 'copy',
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:params.report_outdir, publish_id:report_outdir) }
    
    if ( params.enable_conda ) {
        if (file("${params.dependency_prefix}/panaroo").exists()) {
            conda "${params.dependency_prefix}/panaroo"
        } else {
            conda "${moduleDir}/environment.yml"
        }
    } else {
        conda null
    }
    

    cache 'lenient'
    scratch true
    
    
    input:
    val(gffs)

    output:
    path('pangenome_statistics.txt'), emit: pangenome_summary
    path('gene_presence_absence_roary.csv'), emit: pangenome_csv
    path('gene_presence_absence.Rtab'), emit: pangenome_rtab
    path('panaroo/*')
    path('version_panaroo.txt'), emit: version
    
    script:
    def gffs_str = gffs.join(' ')
    """
    mkdir results
    panaroo -i $gffs_str -o panaroo --clean-mode strict --remove-invalid-genes
    csvtk add-header -t -T -n 'Genes,Range,Total' panaroo/summary_statistics.txt > pangenome_statistics.txt
    cp panaroo/gene_presence_absence_roary.csv .
    cp panaroo/gene_presence_absence.Rtab gene_presence_absence.Rtab
    echo -e panaroo'\t'\$CONDA_PREFIX'\t'\$(panaroo --version)'\t'${params.panaroo_ref} | csvtk add-header -t -n 'tool,conda_env,version,reference' > version_panaroo.txt
    """
}
