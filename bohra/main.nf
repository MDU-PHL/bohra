#!/usr/bin/env nextflow
nextflow.preview.dsl=2
version = '1.0'

println "Ready"
println "Welcome to bohra - microbial genomics pipeline!"

println "Setting up some variables for pipeline."
println "The pipeline will be run in $launchDir"
println "Job directory is ${params.outdir}"
// params.prefill_path = '/home/seq/MDU/QC'
// SAMPLE = config['isolates'].split()
// SNIPPY_SINGULARITY = config['snippy_singularity']
// # ASSEMBLER_SINGULARITY = config['assembler_singularity']
// ABRITAMR_SINGULARITY = config['abritamr_singularity']
// params.min_aln = 0
// params.min_cov = 0
// params.min_qscore = 30
// params.reference_fasta = file('Lm_Cluster1_J1-108.fa' )
// params.outdir = 'TEST_NEXTFLOW'
// params.publish_dir_mode = 'copy'
// // threads
// params.snippy_threads = 4
// params.kraken_threads = 16
// params.assembler_threads = 8
// params.prokka_threads = 8
// params.iqtree_cpus = 20
// // mode
// params.mode = 'sa'
// // template_dir
// params.template_dir = file("${projectDir}/templates")
// println params.template_dir
// // tools to use
// params.run_kraken = true
// params.kraken2_db = "/home/linuxbrew/db/kraken2/pluspf"
// params.assembler = "shovill"
// params.mask_string = ""
params.template_dir = file("${projectDir}/templates")
// println 
def samples = []
input_file = file(params.isolates) // need to make this an input file 
    reader =   input_file.newReader()
    reader.eachLine { line ->
        // print(line)
    samples << line
    }
// println samples
def contigs = [:]
// open the distribution table

contigs_file = file(params.contigs_file) // need to make this an input file 
contigs_reader =   contigs_file.newReader()
// read the file
contigs_reader.eachLine { line ->
    if( line.split('\t')[0] != 'MDUID') {
        // println line.split('\t')
        contigs[line.split('\t')[0]] = line.split('\t')[1]
        }
    }



reads = Channel.fromFilePairs("${params.outdir}/*/*_R{1,2}*.f*q.gz")
                .filter { sample, files -> samples.contains(files[0].getParent().getName())}
                .map { sample, files -> tuple([id: sample, single_end:false, contigs: contigs[sample]], files)}


// println reads.view()

workflow {
    
    include { READ_ANALYSIS;RUN_KRAKEN } from './workflows/common'
    include { PREVIEW_NEWICK } from './workflows/preview'
    include { COLLATE_KRAKEN;COLLATE_SEQS;WRITE_HTML } from './workflows/collation'
    // include { RUN_SNIPPY } from './workflows/snps'
    // include { RUN_CORE } from './workflows/core'
    include { RUN_ASSEMBLE;CONCAT_MLST;CONCAT_RESISTOMES;COLLATE_ASM_PROKKA;CONCAT_ASM;RUN_SNIPPY;RUN_CORE;CONCAT_VIRULOMES } from './workflows/snps'
    

    READ_ANALYSIS ( reads )
    if (params.mode == 'preview') {
        PREVIEW_NEWICK ( READ_ANALYSIS.out.skch.map { cfg, sketch -> sketch }.collect() )
        COLLATE_SEQS ( READ_ANALYSIS.out.stats.map { cfg, stats -> stats }.collect() )
        results = PREVIEW_NEWICK.out.nwk.concat( COLLATE_SEQS.out.collated_seqdata )
    } else if (params.mode == 'default'){
        RUN_SNIPPY ( reads )
        // println RUN_SNIPPY.out.aln.map { cfg, aln -> cfg.id }.collect().view()
        RUN_CORE ( RUN_SNIPPY.out.aln.map { cfg, aln -> cfg.id }.collect() )
        RUN_ASSEMBLE ( reads )
        // println RUN_ASSEMBLE.out.assembly_stats.map { cfg, assembly_stats -> assembly_stats }.collect().map { files -> tuple("assembly", files)}.view()
        // println Channel.value("assembly").view()
        CONCAT_MLST ( RUN_ASSEMBLE.out.mlst.map { cfg, mlst -> mlst }.collect().map { files -> tuple("mlst", files)} )
        CONCAT_RESISTOMES ( RUN_ASSEMBLE.out.resistome.map { cfg, resistome -> resistome }.collect().map { files -> tuple("resistome", files)} )
        CONCAT_VIRULOMES ( RUN_ASSEMBLE.out.virulome.map { cfg, resistome -> resistome }.collect().map { files -> tuple("virulome", files)} )
        // combined asm and prokka stats
        APS = RUN_ASSEMBLE.out.prokka_txt.join( RUN_ASSEMBLE.out.assembly_stats )
        COLLATE_ASM_PROKKA ( APS )
        CONCAT_ASM ( COLLATE_ASM_PROKKA.out.collated_asm.map { cfg, asm -> asm }.collect().map { files -> tuple("assembly", files)} )
        results = CONCAT_ASM.out.collated_assembly.concat( RUN_CORE.out.core_stats, CONCAT_VIRULOMES.out.collated_virulomes, CONCAT_RESISTOMES.out.collated_resistomes, CONCAT_MLST.out.collated_mlst, RUN_CORE.out.newick )

    }

    if ( params.run_kraken ) {
            RUN_KRAKEN ( reads )
            // println RUN_KRAKEN.out.species.map { cfg, species -> species }.collect().view()
            COLLATE_KRAKEN ( RUN_KRAKEN.out.species.map { cfg, species -> species }.collect() )
            results = results.concat( COLLATE_KRAKEN.out.collated_species )
        }
    println results.collect().view()
    // WRITE_HTML ( results.collect() )
    
    
    
                                
                
}