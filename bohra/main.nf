#!/usr/bin/env nextflow
nextflow.preview.dsl=2
version = '1.0'

println "Ready"
println "Welcome to bohra - microbial genomics pipeline!"

println "Setting up some variables for pipeline."
println "The pipeline will be run in $launchDir"


params.prefill_path = '/home/seq/MDU/QC'
// SAMPLE = config['isolates'].split()
// SNIPPY_SINGULARITY = config['snippy_singularity']
// # ASSEMBLER_SINGULARITY = config['assembler_singularity']
// ABRITAMR_SINGULARITY = config['abritamr_singularity']
params.min_aln = 0
params.min_cov = 0
params.min_qscore = 30
params.reference_fasta = file('Lm_Cluster1_J1-108.fa' )
params.outdir = 'TEST_NEXTFLOW'
params.publish_dir_mode = 'copy'
// threads
params.snippy_threads = 4
params.kraken_threads = 16
params.assembler_threads = 8
params.prokka_threads = 8
params.iqtree_cpus = 20
// mode
params.mode = 'sa'
// template_dir
params.template_dir = file("${projectDir}/templates")
println params.template_dir
// tools to use
params.run_kraken = true
params.kraken2_db = "/home/linuxbrew/db/kraken2/pluspf"
params.assembler = "shovill"
params.mask_string = ""

reads = Channel.fromFilePairs("${params.outdir}/*/*_R{1,2}.fq.gz")
                                                                .map { sample, files -> tuple([id: sample, single_end:false], files)}
// println reads.view()


workflow {
    
    include { PREVIEW_ANALYSIS;PREVIEW_NEWICK;RUN_KRAKEN } from './workflows/preview'
    include { COLLATE_KRAKEN;COLLATE_SEQS;WRITE_HTML } from './workflows/collation'
    include { RUN_SNIPPY } from './workflows/snps'
    include { RUN_CORE } from './workflows/core'
    include { RUN_ASSEMBLE;CONCAT_MLST;CONCAT_RESISTOMES } from './workflows/assemble_typing'
    

    PREVIEW_ANALYSIS ( reads )
    PREVIEW_NEWICK ( PREVIEW_ANALYSIS.out.skch.map { cfg, sketch -> sketch }.collect() )
    COLLATE_SEQS ( PREVIEW_ANALYSIS.out.stats.map { cfg, stats -> stats }.collect() )
    preview_results = PREVIEW_NEWICK.out.nwk.concat( COLLATE_SEQS.out.collated_seqdata )
    if ( params.run_kraken ) {
        RUN_KRAKEN ( reads )
        // println RUN_KRAKEN.out.species.map { cfg, species -> species }.collect().view()
        COLLATE_KRAKEN ( RUN_KRAKEN.out.species.map { cfg, species -> species }.collect() )
        preview_results = preview_results.concat( COLLATE_KRAKEN.out.collated_species )
    }
    
    if (params.mode == 'sa'){
        RUN_SNIPPY ( reads )
        // println RUN_SNIPPY.out.aln.map { cfg, aln -> cfg.id }.collect().view()
        RUN_CORE ( RUN_SNIPPY.out.aln.map { cfg, aln -> cfg.id }.collect() )
        RUN_ASSEMBLE ( reads )
        println RUN_ASSEMBLE.out.assembly_stats.map { cfg, assembly_stats -> assembly_stats }.collect().map { files -> tuple("assembly", files)}.view()
        // println Channel.value("assembly").view()
        CONCAT_MLST ( RUN_ASSEMBLE.out.mlst.map { cfg, mlst -> mlst }.collect().map { files -> tuple("mlst", files)} )
        CONCAT_RESISTOMES ( RUN_ASSEMBLE.out.resistome.map { cfg, resistome -> resistome }.collect().map { files -> tuple("resistome", files)} )
        // COLLATE_MLST ( RUN_ASSEMBLE.out.mlst )
        // COLLA
    }

    // WRITE_HTML ( preview_results.collect() )
    
    
    
                                
                
}