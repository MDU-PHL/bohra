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
    // publishDir "abritamr",
    //     mode: 'link',
    //     saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/abritamr" : 'bioconda::abritamr') : null) 
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}/bohra-abritamr").exists()) {
            conda "${params.conda_path}/bohra-abritamr"
        } else {
            conda 'environment.yml'
        }
    } else {
        conda null
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
    echo -e abritamr'\t'\$CONDA_PREFIX'\t'\$(abritamr -v) | csvtk add-header -t -n 'tool,conda_env,version' > version_abritamr.txt
    """
}

process ABRITAMR_GENERAL {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"${meta.id}", publish_id:meta.id) }
    // publishDir "abritamr",
    //     mode: 'link',
    //     saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/abritamr" : 'bioconda::abritamr') : null) 
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-abritamr"
        } else {
            conda 'bioconda::bioconda::abritamr'
        }
    } else {
        conda null
    }

    scratch true
    input:
    tuple val(meta), path(summary_matches), path(summary_partials)

    output:
    tuple val(meta), path("reportable_amr_matches.txt"), emit: reportable

    """
    echo "ISOLATE,SPECIES_EXP,SPECIES_OBS,TEST_QC\n${meta.id},${meta.species},${meta.species},PASS" > qc.tmp.csv
    abritamr report -q qc.tmp.csv -r bohra -m $summary_matches -p $summary_partials --sop general --sop_name reportable
    csvtk xlsx2csv -n reportable bohra_reportable.xlsx | csvtk rename -f 'MDU sample ID,Resistance genes (alleles) detected,Resistance genes (alleles) det (non-rpt)' -n 'Isolate,Reportable AMR mechanisms,Other AMR mechanisms' | csvtk cut -f 'Isolate,Reportable AMR mechanisms,Other AMR mechanisms' |  csvtk csv2tab > reportable_amr_matches.txt
    """
    
}



process ABRITAMR_INFER {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"${meta.id}", publish_id:meta.id) }
    // publishDir "abritamr",
    //     mode: 'link',
    //     saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id, publish_id:meta.id) }
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/abritamr" : 'bioconda::abritamr') : null) 
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-abritamr"
        } else {
            conda 'bioconda::bioconda::abritamr'
        }
    } else {
        conda null
    }
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

// process COLLATE_ABRITAMR {
    
//     label 'process_low'
        
//     cache 'lenient'
//     input:
//     val(abritamr_matches_isolates) 
    
//     output:
//     path 'summary_matches.txt', emit: abritamr_matches
    
    
//     script:
    
//     """
//     $module_dir/concat.py summary_matches.txt $abritamr_matches_isolates
//     """
    
// }

// process COLLATE_ABRITAMR_PARTIALS {
    
//     label 'process_low'
      
//     cache 'lenient'
//     input:
//     val(abritamr_partials_isolates)

//     output:
//     path 'summary_partials.txt', emit: abritamr_partials
    
//     script:
    
//     """
//     $module_dir/concat.py summary_partials.txt $abritamr_partials_isolates
//     """
    
// }
