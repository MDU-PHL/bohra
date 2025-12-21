// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

module_dir = moduleDir + "/bin"

params.options = [:]
def options    = initOptions(params.options)


process ABRITAMR {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"${meta.id}", publish_id:meta.id) }
    
    if ( params.enable_conda ) {
        if (file("${params.dependency_prefix}/abritamr").exists()) {
            conda "${params.dependency_prefix}/abritamr"
        }  
        }
    

    cache 'lenient' 
    scratch true
    // afterScript "rm -fr /tmp/\$USER/*"
    
    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("summary_virulence.txt"), emit: abritamr_virulence
    tuple val(meta),path('summary_matches.txt'), emit: abritamr_matches
    tuple val(meta),path('summary_partials.txt'), emit: abritamr_partials
    tuple val(meta),path('amrfinder.out'), emit: amrfinder_out
    tuple val(meta), path('version_abritamr.txt'), emit: version
    script:
    
    
    """
    sp=\$($module_dir/extract_species.py $contigs)
    abritamr run -c $contigs -px ${meta.id} -j $task.cpus \$sp
    cp ${meta.id}/* .
    $module_dir/add_species.py '${meta.species}' summary_matches.txt
    echo -e abritamr'\t'\$CONDA_PREFIX'\t'\$(abritamr -v)'\t'${params.abritamr_ref} | csvtk add-header -t -n 'tool,conda_env,version,reference' > version_abritamr.txt
    """
}

process ABRITAMR_GENERAL {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"${meta.id}", publish_id:meta.id) }
    if (file("${params.dependency_prefix}/abritamr").exists()) {
            conda "${params.dependency_prefix}/abritamr"
        } 
    scratch true
    input:
    tuple val(meta), path(summary_matches), path(summary_partials)

    output:
    tuple val(meta), path("reportable_amr_matches.txt"), emit: reportable

    """
    echo "ISOLATE,SPECIES_EXP,SPECIES_OBS,TEST_QC\n${meta.id},${meta.species},${meta.species},PASS" > qc.tmp.csv
    abritamr report -q qc.tmp.csv -r bohra -m $summary_matches -p $summary_partials --sop general --sop_name reportable
    csvtk xlsx2csv -n reportable bohra_reportable.xlsx | csvtk rename -f 'MDU sample ID,Resistance genes (alleles) detected,Resistance genes (alleles) det (non-rpt),Species_obs' -n 'Isolate,Reportable AMR mechanisms,Other AMR mechanisms,Species' | csvtk cut -f 'Isolate,Reportable AMR mechanisms,Other AMR mechanisms,Species' |  csvtk csv2tab > reportable_amr_matches_tmp.txt
    $module_dir/separate_mechanisms.py reportable_amr_matches_tmp.txt
    """
    
}



process ABRITAMR_INFER {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"${meta.id}", publish_id:meta.id) }
    if (file("${params.dependency_prefix}/abritamr").exists()) {
            conda "${params.dependency_prefix}/abritamr"
        } 
    errorStrategy 'ignore'
    scratch true
    input:
    tuple val(meta), path(summary_matches)
    

    output:
    tuple val(meta), path("bohra_inferred.csv"), emit: inferred
    
    script:
    
    """
    echo "ISOLATE,SPECIES_EXP,SPECIES_OBS,TEST_QC\n${meta.id},${meta.species},${meta.species},PASS" > qc.tmp.csv
    abritamr report -q qc.tmp.csv -r bohra -m $summary_matches --sop plus --sop_name inferred
    csvtk xlsx2csv -n bohra_inferred.xlsx | csvtk csv2tab > bohra_inferred.csv
    """
    
}



process COMBINE_AMR {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    
    cache 'lenient' 
    scratch true

    input:
    tuple val(meta), path(matches), path(partials), path(amrfinder_out), path(plasmids)

    output:
    tuple val(meta), path("resistome.txt"), emit: resistome
    tuple val(meta), path("plasmid.json"), emit: plasmid
    
    script:
    """
    $module_dir/combine.py ${meta.id} $matches $partials $amrfinder_out $plasmids
    """

}

