#!/usr/bin/env nextflow
nextflow.enable.dsl=2
version = '1.0'

println "Ready"
println "Running bohra - microbial genomics pipeline!"

println "Setting up some variables for pipeline."
println "The pipeline will be run in $launchDir"
println "Job directory is ${params.outdir}"
println "The conda path is ${params.conda_path}"
e = file("${params.conda_path}").exists()
println "The conda path exists : $e"

println "Conda has been enabled : ${params.enable_conda}"
println "The blast_db is : ${params.blast_db}"
// set some parameters
params.species_options = ['Neisseria', 'Acinetobacter_baumannii', "Campylobacter", "Enterococcus_faecalis", "Enterococcus_faecium", "Escherichia", "Klebsiella", "Salmonella", "Staphylococcus_aureus", "Staphylococcus_pseudintermedius", "Streptococcus_agalactiae", "Streptococcus_pneumoniae", "Streptococcus_pyogenes", "Vibrio_cholerae"]
params.template_dir = file("${projectDir}/templates")

def samples = []
input_file = file(params.isolates) // need to make this an input file 
    reader =   input_file.newReader()
    reader.eachLine { line ->
    samples << line
    }
def contigs = [:]
// open the distribution table
if (params.contigs_file != 'no_contigs'){
contigs_file = file(params.contigs_file) // need to make this an input file 
contigs_reader =   contigs_file.newReader()
// read the file
contigs_reader.eachLine { line ->
    if( line.split('\t')[0] != 'MDUID') {
        contigs[line.split('\t')[0]] = line.split('\t')[1]
        }
    }
} 

// println contigs

reference = Channel.fromPath( "${params.reference}")
reads = Channel.fromFilePairs(["${params.outdir}/*/R{1,2}*.f*q.gz","${params.outdir}/*/*_{1,2}.f*q.gz"])
                .filter { sample, files -> samples.contains(files[0].getParent().getName())}
                .map { sample, files -> tuple([id: files[0].getParent().getName(), single_end:false, contigs: contigs[files[0].getParent().getName()]], files)}





include { READ_ANALYSIS;RUN_KRAKEN } from './workflows/common'
include { PREVIEW_NEWICK } from './workflows/preview'
include { COLLATE_KRAKEN;COLLATE_SEQS;WRITE_HTML } from './workflows/collation'
// include { RUN_SNIPPY } from './workflows/snps'
include { RUN_PANAROO } from './workflows/pangenome'
include { RUN_ASSEMBLE;CONCAT_STATS;CONCAT_MLST;CONCAT_RESISTOMES;COLLATE_ASM_PROKKA;CONCAT_ASM;RUN_SNIPPY;RUN_CORE;RUN_IQTREE;CONCAT_VIRULOMES;CONCAT_CORE_STATS;CONCAT_PLASMID;RUN_GUBBINS } from './workflows/default'
    
workflow {
    
    

    READ_ANALYSIS ( reads,reference )
    if (params.mode == 'preview') {
        PREVIEW_NEWICK ( READ_ANALYSIS.out.skch.map { cfg, sketch -> sketch }.collect() )
        COLLATE_SEQS ( READ_ANALYSIS.out.stats.map { cfg, stats -> stats }.collect() )
        results = PREVIEW_NEWICK.out.nwk.concat( COLLATE_SEQS.out.collated_seqdata )
    } else if (params.mode != 'preview'){
        RUN_SNIPPY ( reads.combine( reference ) )
        RUN_CORE ( RUN_SNIPPY.out.aln.map { cfg, aln -> aln.getParent() }.collect(), reference )
        core_aln =  RUN_CORE.out.core_aln
        if ( params.gubbins ){
            RUN_GUBBINS( RUN_CORE.out.core_full_aln )
            core_aln = RUN_GUBBINS.out.core_aln
        }
        if (params.run_iqtree ){
            RUN_IQTREE ( core_aln, RUN_CORE.out.core_full_aln)
            tree = RUN_IQTREE.out.newick
        } else {
            tree = Channel.empty().ifEmpty('EmptyFile')
        }
        RUN_ASSEMBLE ( reads )
        if (params.mode == 'pluspan') {
            RUN_PANAROO( RUN_ASSEMBLE.out.gff.map { cfg, gff -> gff }.collect() )
            svg = RUN_PANAROO.out.svg
        } else {
            svg = Channel.empty().ifEmpty('EmptyFile')
        }
        CONCAT_MLST ( RUN_ASSEMBLE.out.mlst.map { cfg, mlst -> mlst }.collect().map { files -> tuple("mlst", files)} )
        CONCAT_STATS ( READ_ANALYSIS.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple("seqdata", files)} )
        CONCAT_RESISTOMES ( RUN_ASSEMBLE.out.resistome.map { cfg, resistome -> resistome }.collect().map { files -> tuple("resistome", files)} )
        CONCAT_VIRULOMES ( RUN_ASSEMBLE.out.virulome.map { cfg, resistome -> resistome }.collect().map { files -> tuple("virulome", files)} )
        // combined asm and prokka stats
        APS = RUN_ASSEMBLE.out.prokka_txt.join( RUN_ASSEMBLE.out.assembly_stats )
        COLLATE_ASM_PROKKA ( APS )
        CONCAT_CORE_STATS ( RUN_SNIPPY.out.qual.map { cfg, core_stats -> core_stats }.collect().map { files -> tuple("core_genome", files)} )
        CONCAT_ASM ( COLLATE_ASM_PROKKA.out.collated_asm.map { cfg, asm -> asm }.collect().map { files -> tuple("assembly", files)} )
        CONCAT_PLASMID ( RUN_ASSEMBLE.out.plasmid.map { cfg, plasmid -> plasmid }.collect().map { files -> tuple("plasmid", files)} )
        results = CONCAT_ASM.out.collated_assembly.concat(svg, tree, CONCAT_PLASMID.out.collated_plasmids, CONCAT_CORE_STATS.out.collated_core, CONCAT_VIRULOMES.out.collated_virulomes, CONCAT_RESISTOMES.out.collated_resistomes, CONCAT_MLST.out.collated_mlst )

    }

    if ( params.run_kraken ) {
            kraken = Channel.fromPath( "${params.kraken2_db}")
            RUN_KRAKEN ( reads.combine(kraken) )
            COLLATE_KRAKEN ( RUN_KRAKEN.out.species.map { cfg, species -> species }.collect() )
            results = results.concat( COLLATE_KRAKEN.out.collated_species )
        }
    WRITE_HTML ( results.collect() ) 
}