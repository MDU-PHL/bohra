// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process ASSEMBLER_PE {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"${meta.id}", publish_id:meta.id) }
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}/bohra-shovill").exists()) {
            conda "${params.conda_path}/bohra-shovill"
        } else {
            conda 'environment.yml'
        }
    } else {
        conda null
    }

    cache 'lenient'
    scratch true
    
    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('contigs.fa') , emit: contigs
    tuple val(meta), path('version_assembler.txt') , emit: version

    script:
    
    if ( meta.asm != "not_supplied")
        """
        if [ -e $meta.asm ]
        then
            cp $meta.asm contigs.fa
            echo -e Assembly file supplied'\t'Not Applicable'\t'${meta.asm} | csvtk add-header -t -n 'tool,conda_env,version' > version_assembler.txt
        else
            echo "No assembly file found" >> assembly.log
        fi
        """
    else if ( meta.asm == "not_supplied" && params.assembler == "shovill" ) {
        """
        shovill --R1 ${reads[0]} --R2 ${reads[1]} --outdir current --cpus $task.cpus --ram 16 
        cp current/contigs.fa contigs.fa
        version=\$(shovill --version)
        echo -e shovill'\t'\$CONDA_PREFIX'\t'\$(shovill -v) | csvtk add-header -t -n 'tool,conda_env,version' > version_assembler.txt
        """
    } else if (meta.asm == "not_supplied" && params.assembler == "spades"  ) {
        """
        tmp_dir=\$(mktemp -d)
        spades.py -1 ${reads[0]} -2 ${reads[1]} -o current -t $task.cpus $options.args2 --tmp-dir \$tmp_dir
        cp current/contigs.fasta contigs.fa
        rm -rf \$tmp_dir
        echo -e spades'\t'\$CONDA_PREFIX'\t'\$(spades.py -v) | csvtk add-header -t -n 'tool,conda_env,version' > version_spades.txt
        """
    } else if (meta.asm == "not_supplied" && params.assembler == "skesa"  ) {
        """
        skesa --fastq ${reads[0]},${reads[1]} --cores $task.cpus --vector_percent 1.0 \
        --contigs_out contigs.fa
        echo -e skesa'\t'\$CONDA_PREFIX'\t'\$(skesa -v) | csvtk add-header -t -n 'tool,conda_env,version' > version_assembler.txt
        """
    } else {
        """
        echo "No assembly performed" >> assembly.log
        """
    }
    
}
