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
// println samples
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


include { READ_ANALYSIS } from './workflows/read_assessment'
include { RUN_KRAKEN } from './workflows/species'
include { PREVIEW_NEWICK } from './workflows/preview'
include { COLLATE_KRAKEN;COLLATE_SEQS;WRITE_HTML } from './workflows/collation'
include { RUN_SNPS } from './workflows/snps'
include { RUN_PANAROO } from './workflows/pangenome'
include { RUN_ASSEMBLE } from './workflows/assemble'
include { BASIC_TYPING;SEROTYPES;CONCAT_TYPER;CONCAT_RESISTOMES;CONCAT_MLST;CONCAT_VIRULOMES;CONCAT_PLASMID } from './workflows/typing'
include { PREVIEW_VERSIONS;FULL_VERSIONS;DEFAULT_VERSIONS;AMR_TYPING_VERSIONS;SNPS_VERSIONS;ASSEMBLE_VERSIONS } from './workflows/versions'

workflow {
    
    

    READ_ANALYSIS ( reads )
    results = COLLATE_SEQS ( READ_ANALYSIS.out.stats.map { cfg, stats -> stats }.collect() )
    if ( params.run_kraken ) {
            kraken = Channel.fromPath( "${params.kraken2_db}")
            RUN_KRAKEN ( reads.combine(kraken) )
            COLLATE_KRAKEN ( RUN_KRAKEN.out.species.map { cfg, species -> species }.collect() )
            results = results.concat( COLLATE_KRAKEN.out.collated_species )
            species = RUN_KRAKEN.out.species_obs
            
        } else {
            species = Channel.empty().ifEmpty('EmptyFile')
        }
    results = results.concat(species)
    if (params.mode == 'preview') {
        PREVIEW_NEWICK ( reads )
        results = results.concat( PREVIEW_NEWICK.out.nwk )
        PREVIEW_VERSIONS()
        results = results.concat( PREVIEW_VERSIONS.out.versions )
        WRITE_HTML ( results.collect() )
    }
    if ( params.mode == 'snps' || params.mode == 'phylogeny') {
        RUN_SNPS ( reads,reference )
        results = results.concat( RUN_SNPS.out.core_stats, RUN_SNPS.out.tree )
        SNPS_VERSIONS( )
        results = results.concat( SNPS_VERSIONS.out.versions)
        WRITE_HTML ( results.collect() )
    } 
    if ( params.mode == 'assemble' || params.mode == 'amr_typing' ){
            RUN_ASSEMBLE ( reads )
            results = results.concat( RUN_ASSEMBLE.out.assembly_stats )
            if ( params.mode == 'assemble' ){
                ASSEMBLE_VERSIONS( )
                results = results.concat( ASSEMBLE_VERSIONS.out.versions )
            }
            if ( params.mode == 'amr_typing'){
                    BASIC_TYPING ( RUN_ASSEMBLE.out.contigs )
                    results = results.concat( BASIC_TYPING.out.mlst, BASIC_TYPING.out.resistome, BASIC_TYPING.out.virulome, BASIC_TYPING.out.plasmid )
                    if ( params.run_kraken ){
                        typing_input = RUN_ASSEMBLE.out.contigs.join( species )
                        SEROTYPES ( typing_input )
                        CONCAT_TYPER ( SEROTYPES.out.typers)
                        results = results.concat(CONCAT_TYPER.out.collated_typers)
                        AMR_TYPING_VERSIONS()
                        results = results.concat( AMR_TYPING_VERSIONS.out.versions )
                    }
                    WRITE_HTML ( results.collect() )
                } else {
                    WRITE_HTML ( results.collect() )
                }
    } 
    if ( params.mode == 'default' || params.mode == 'full') {
        RUN_SNPS ( reads,reference )
        results = results.concat( RUN_SNPS.out.core_stats, RUN_SNPS.out.tree )
        RUN_ASSEMBLE ( reads )
        results = results.concat( RUN_ASSEMBLE.out.assembly_stats )
        BASIC_TYPING ( RUN_ASSEMBLE.out.contigs )
        results = results.concat( BASIC_TYPING.out.mlst, BASIC_TYPING.out.resistome, BASIC_TYPING.out.virulome, BASIC_TYPING.out.plasmid )
        if ( params.mode == 'default' ) {
            DEFAULT_VERSIONS()
            results = results.concat( DEFAULT_VERSIONS.out.versions )
        } else if (params.mode == 'full') {
            FULL_VERSIONS()
            results = results.concat( FULL_VERSIONS.out.versions )
        }
        if ( params.run_kraken ){
            typing_input = RUN_ASSEMBLE.out.contigs.join( species )
            SEROTYPES ( typing_input )
            CONCAT_TYPER ( SEROTYPES.out.typers)
            results = results.concat(CONCAT_TYPER.out.collated_typers)
        }
        if (params.mode == 'full') {
                    RUN_PANAROO( RUN_ASSEMBLE.out.gff.map { cfg, gff -> gff }.collect() )
                    svg = RUN_PANAROO.out.svg
                    
                } else {
                    svg = Channel.empty().ifEmpty('EmptyFile')
                }
        results = results.concat( svg )
        WRITE_HTML ( results.collect() )
    }

    
     
}