// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'
module_dir = moduleDir + "/bin"
params.options = [:]
def options    = initOptions(params.options)



process EXTRACT_SPECIES {
    tag "$meta.id"
    label 'process_low'

    input:
        tuple val(meta), path(kraken2)
    output:
        tuple val(meta), env(species), emit: species_obs

    script:
    """
    species=\$($module_dir/extract_species.py  $kraken2)
    """

}


process VERSION_KMC {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-kmc"
        } else {
            conda 'kmc csvtk'
        }
    } else {
        conda null
    }
    output:
    path 'version_kmc.txt', emit: version

    script:

    """
    echo -e kmc'\t'\$CONDA_PREFIX'\t'\$(kmc -V | grep 'K-Mer Counter') | csvtk add-header -t -n 'tool,conda_env,version' > version_kmc.txt
    """
}


process VERSION_EMMTYPER {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-emmtyper"
        } else {
            conda 'emmtyper csvtk'
        }
        
    } else {
        conda null
    }
    output:
    path "version_emmtyper.txt", emit: version

    script:
    """
    echo -e emmtyper'\t'\$CONDA_PREFIX'\t'\$(emmtyper --version) | csvtk add-header -t -n 'tool,conda_env,version' > version_emmtyper.txt
    """
}

process VERSION_ECTYPER {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-ectyper"
        } else {
            conda 'ectyper csvtk'
        }
        
    } else {
        conda null
    }
    output:
    path "version_ectyper.txt", emit: version

    script:
    """
    echo -e ectyper'\t'\$CONDA_PREFIX'\t'\$(ectyper --version) | csvtk add-header -t -n 'tool,conda_env,version' > version_ectyper.txt
    """
}

process VERSION_KLEBORATE {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-kleborate"
        } else {
            conda 'kleborate csvtk'
        }
        
    } else {
        conda null
    }
    output:
    path "version_kleborate.txt", emit: version

    script:
    """
    echo -e kleborate'\t'\$CONDA_PREFIX'\t'\$(kleborate --version) | csvtk add-header -t -n 'tool,conda_env,version' > version_kleborate.txt
    """
}

process VERSION_NGMASTER {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-ngmaster"
        } else {
            conda 'ngmaster csvtk'
        }
       
    } else {
        conda null
    }
    output:
    path "version_ngmaster.txt", emit: version

    script:
    """
    echo -e ngmaster'\t'\$CONDA_PREFIX'\t'\$(ngmaster --version) | csvtk add-header -t -n 'tool,conda_env,version' > version_ngmaster.txt
    """
}

process VERSION_MENINGOTYPE {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-meningotype"
        } else {
            conda 'meningotype csvtk'
        }
    } else {
        conda null
    }
    output:
    path "version_meningotype.txt", emit: version

    script:
    """
    echo -e meningotype'\t'\$CONDA_PREFIX'\t'\$(meningotype --version) | csvtk add-header -t -n 'tool,conda_env,version' > version_meningotype.txt
    """
}

process VERSION_LISSERO {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-lissero"
        } else {
            conda 'lissero csvtk'
        }
        
    } else {
        conda null
    }

    output:
    path "version_lissero.txt", emit: version

    script:
    """
    echo -e lissero'\t'\$CONDA_PREFIX'\t'\$(lissero --version) | csvtk add-header -t -n 'tool,conda_env,version' > version_lissero.txt
    """
}



process VERSION_STYPE {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-stype"
        } 
        // will need to release stype to conda added in ignore strategy in case people don't use init - at least whole pipeline won't fall down
    } else {
        conda null
    }

    output:
    path "version_stype.txt", emit: version

    script:
    """
    echo -e stype'\t'\$CONDA_PREFIX'\t'\$(stype -v) | csvtk add-header -t -n 'tool,conda_env,version' > version_stype.txt
    """
}



process VERSION_MOBSUITE {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-mob_suite"
        } else {
            conda 'bioconda::mob_suite=3.0.2 csvtk'
        }
    } else {
        conda null
    }

    output:
    path "version_mobsuite.txt", emit: version

    script:
    """
    echo -e mobsuite'\t'\$CONDA_PREFIX'\t'\$(mob_recon -V) | csvtk add-header -t -n 'tool,conda_env,version' > version_mobsuite.txt
    """
}

process VERSION_QUICKTREE {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-quicktree"
        } else {
            conda 'bioconda::bioconda::quicktree=2.5 newick_utils csvtk'
        }
    } else {
        conda null
    }

    output:
    path "version_quicktree.txt", emit: version

    script:
    """
    echo -e quicktree'\t'\$CONDA_PREFIX'\t'\$(quicktree -v) | csvtk add-header -t -n 'tool,conda_env,version' > version_quicktree.txt
    """
}

process VERSION_MASH {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-mash"
        } else {
            conda 'mash csvtk'
        }
    } else {
        conda null
    }

    output:
    path "version_mash.txt", emit: version

    script:
    """
    echo -e mash'\t'\$CONDA_PREFIX'\t'\$(mash --version) | csvtk add-header -t -n 'tool,conda_env,version' > version_mash.txt
    """
}

process VERSION_PROKKA {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-prokka"
        } else {
            conda 'bioconda::prokka csvtk'
        }
    } else {
        conda null
    }

    output:
    path "version_prokka.txt", emit: version

    script:
    """
    echo -e prokka'\t'\$CONDA_PREFIX'\t'\$(prokka -v) | csvtk add-header -t -n 'tool,conda_env,version' > version_prokka.txt
    """
}


process VERSION_MLST {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-mlst"
        } else {
            conda 'bioconda::mlst=2.19.0 csvtk'
        }
    } else {
        conda null
    }

    output:
    path "version_mlst.txt", emit: version

    script:
    """
    echo -e mlst'\t'\$CONDA_PREFIX'\t'\$(mlst -v) | csvtk add-header -t -n 'tool,conda_env,version' > version_mlst.txt
    """
}

process VERSION_ABRITAMR {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-abritamr"
        } else {
            conda 'bioconda::bioconda::abritamr csvtk'
        }
    } else {
        conda null
    }

    output:
    path "version_abritamr.txt", emit: version

    script:
    """
    echo -e abritamr'\t'\$CONDA_PREFIX'\t'\$(abritamr -v) | csvtk add-header -t -n 'tool,conda_env,version' > version_abritamr.txt
    """
}

process VERSION_SPADES {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-spades"
        } else {
            conda 'bioconda::spades=3.15.2 csvtk'
        }
    } else {
        conda null
    }

    output:
    path "version_assembler.txt", emit: version

    script:
    if( params.contigs_file == 'no_contigs')
    """
    echo -e spades'\t'\$CONDA_PREFIX'\t'\$(spades.py -v) | csvtk add-header -t -n 'tool,conda_env,version' > version_spades.txt
    """
    else
    """
    echo -e Assembly file supplied'\t'Not Applicable'\t'${params.contigs_file} | csvtk add-header -t -n 'tool,conda_env,version' > version_assembler.txt
    """
}

process VERSION_SKESA {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-skesa"
        } else {
            conda 'bioconda::skesa=2.4 csvtk'
        }
    } else {
        conda null
    }

    output:
    path "version_assembler.txt", emit: version

    script:
    if( params.contigs_file == 'no_contigs')
    """
    echo -e skesa'\t'\$CONDA_PREFIX'\t'\$(skesa -v) | csvtk add-header -t -n 'tool,conda_env,version' > version_assembler.txt
    """
    else
    """
    echo -e Assembly file supplied'\t'Not Applicable'\t'${params.contigs_file} | csvtk add-header -t -n 'tool,conda_env,version' > version_assembler.txt
    """
}


process VERSION_SHOVILL {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-shovill"
        } else {
            conda 'bioconda::shovill=1.1.0 csvtk'
        }
    } else {
        conda null
    }
    output:
    path "version_assembler.txt", emit: version

    script: 
    if( params.contigs_file == 'no_contigs' )
        """
        echo -e shovill'\t'\$CONDA_PREFIX'\t'\$(shovill -v) | csvtk add-header -t -n 'tool,conda_env,version' > version_assembler.txt
        """
    else
        """
        echo -e Assembly file supplied'\t'Not Applicable'\t'${params.contigs_file} | csvtk add-header -t -n 'tool,conda_env,version' > version_assembler.txt
        """
}


process VERSION_PANAROO {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-panaroo"
        } else {
            conda 'bioconda::panaroo=1.2.9 csvtk'
        }
    } else {
        conda null
    }
    output:
    path "version_panaroo.txt", emit: version

    script:
    """
    echo -e panaroo'\t'\$CONDA_PREFIX'\t'\$(panaroo --version ) | csvtk add-header -t -n 'tool,conda_env,version' > version_panaroo.txt
    """
}


process VERSION_GUBBINS {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-gubbins"
        } else {
            conda 'gubbins=2.4.1 snp-sites=2.5.1 csvtk'
        }
    } else {
        conda null
    }

    output:
    path "version_gubbins.txt", emit: version

    script:
    """
    echo -e gubbins'\t'\$CONDA_PREFIX'\t'\$(run_gubbins.py --version) | csvtk add-header -t -n 'tool,conda_env,version' > version_gubbins.txt
    """
}



process VERSION_ANY2FASTA {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    output:
    path "version_any2fasta.txt", emit: version

    script:
    """
    echo -e any2fasta'\t'\$CONDA_PREFIX'\t'\$(any2fasta -v) | csvtk add-header -t -n 'tool,conda_env,version' > version_any2fasta.txt
    """
}


process VERSION_SEQKIT {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-seqkit"
        } else {
            conda 'csvtk seqkit=2.1.0'
        }
    } else {
        conda null
    }

    output:
    path "version_seqkit.txt", emit: version

    script:
    """
    echo -e seqkit'\t'\$CONDA_PREFIX'\t'\$(seqkit version) | csvtk add-header -t -n 'tool,conda_env,version' > version_seqkit.txt
    """
}


process VERSION_KRAKEN2 {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-kraken2"
        } else {
            conda 'kraken2=2.1.2 csvtk'
        }
    } else {
        conda null
    }

    output:
    path "version_kraken2.txt", emit: version

    script:
    """
    echo -e kraken2'\t'\$CONDA_PREFIX'\t'\$(kraken2 --version | grep version) | csvtk add-header -t -n 'tool,conda_env,version' > version_kraken2.txt
    """
}

process VERSION_NEXTFLOW {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    output:
    path "version_nextflow.txt", emit: version

    script:
    """
    echo -e Nextflow'\t'\$CONDA_PREFIX'\t'\$(nextflow -v) | csvtk add-header -t -n 'tool,conda_env,version' > version_nextflow.txt
    """
}


process VERSION_SNIPPY {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-snippy"
        } else {
            conda 'bioconda::snippy=4.4.5 csvtk'
        }
    } else {
        conda null
    }


    output:
    path "version_snippy.txt", emit: version

    script:
    """
    echo -e snippy'\t'\$CONDA_PREFIX'\t'\$(snippy -v 2>&1) | csvtk add-header -t -n 'tool,conda_env,version' > version_snippy.txt
    """
}


process VERSION_SNPDISTS {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-snpdists"
        } else {
            conda 'bioconda::snp-dists=0.8.2 bioconda::csvtk'
        }
    } else {
        conda null
    }
    
    output:
    path "version_snpdists.txt", emit: version

    script:
    """
    echo -e snp-dists'\t'\$CONDA_PREFIX'\t'\$(snp-dists -v) | csvtk add-header -t -n 'tool,conda_env,version' > version_snpdists.txt
    """
}



process VERSION_IQTREE {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', publish_id:'report') }
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-iqtree"
        } else {
            conda 'iqtree=2.1.4 snp-sites=2.5.1'
        }
    } else {
        conda null
    }

    output:
    path "version_iqtree.txt", emit: version

    script:
    """
    echo -e iQtree'\t'\$CONDA_PREFIX'\t'\$(iqtree --version | grep version) | csvtk add-header -t -n 'tool,conda_env,version' > version_iqtree.txt
    """
}