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
    
    
    // conda (params.enable_conda ? (file("${params.conda_path}").exists() ? "${params.conda_path}/shovill" : 'shovill=1.1.0') : null) 
    
    if ( params.enable_conda ) {
        if (file("${params.conda_path}").exists()) {
            conda "${params.conda_path}/bohra-shovill"
        } else {
            conda 'bioconda::shovill=1.1.0'
        }
    } else {
        conda null
    }

    cache 'lenient'
    // afterScript "rm -fr /tmp/\$USER/*"
    // scratch true
    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('contigs.fa') , emit: contigs
    tuple val(meta), path('assembly.log') , emit: log

    script:
    
    if ( meta.asm != "no_contigs")
        """
        if [ -e $meta.asm ]
        then
            cp $meta.asm contigs.fa
            echo "Using existing assembly file" >> assembly.log
        else
            echo "No assembly file found" >> assembly.log
        fi
        """
    else if ( meta.asm == "no_contigs" && meta.assembler == "shovill" ) {
        """
        shovill --R1 ${reads[0]} --R2 ${reads[1]} --outdir current --cpus $task.cpus --ram 16 
        cp current/contigs.fa contigs.fa
        version=\$(shovill --version)
        echo "Assembled with shovill version: \$version" >> assembly.log
        cat current/shovill.log >> assembly.log
        """
    } else if (meta.asm == "no_contigs" && meta.assembler == "spades"  ) {
        """
        tmp_dir=\$(mktemp -d)
        spades.py -1 ${reads[0]} -2 ${reads[1]} -o current -t $task.cpus $options.args2 --tmp-dir \$tmp_dir
        cp current/contigs.fasta contigs.fa
        rm -rf \$tmp_dir
        version=\$(spades.py --version)
        echo "Assembled with spades version: \$version" >> assembly.log
        """
    } else if (meta.asm == "no_contigs" && meta.assembler == "skesa"  ) {
        """
        skesa --fastq ${reads[0]},${reads[1]} --cores $task.cpus --vector_percent 1.0 \
        --contigs_out contigs.fa
        version=\$(skesa --version)
        echo "Assembled with skesa version: \$version" >> assembly.log
        """
    } else {
        """
        echo "No assembly performed" >> assembly.log
        """
    }
    
}
