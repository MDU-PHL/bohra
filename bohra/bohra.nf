#!/usr/bin/env nextflow
nextflow.enable.dsl=2
version = '1.0'

println "Ready"
println "Running bohra - microbial genomics pipeline!"

// println "Setting up some variables for pipeline."
// println "The pipeline will be run in $launchDir"
// println "${params.workflows} workflow is being run"
// println "Job directory is ${params.outdir}"
// println "The conda path is ${params.conda_path}"
// e = file("${params.conda_path}").exists()
// println "The conda path exists : $e"

// println "Conda has been enabled : ${params.enable_conda}"
// println "The blast_db is : ${params.blast_db}"
// set some parameters
// params.template_dir = file("${projectDir}/templates")

// get a list of samples
def samples = [:]
input_file = file(params.isolates) // need to make this an input file 
    reader =   input_file.newReader()
    reader.eachLine { line ->
    line = line.split("\t")
        samples[line[0]] = line[1]
    }
// get reads
println "The samples are : $samples"
reads_pe = Channel.fromFilePairs(["${params.outdir}/*/R{1,2}.fastq.gz"])
                .filter { sample, files -> samples.containsKey(files[0].getParent().getName())}
                .map { sample, files -> tuple([id: files[0].getParent().getName(),modules: samples[files[0].getParent().getName()],asm:false, pe_reads:true], files)}
// get assemblies

// println "The reads are : $reads_pe"
// println "The assemblies are : $asm"

// println  reads_pe.view()
// println asm.view()
// println joint_input.view()
// "speciation",
//     "snps",
//     "ska",
//     "tree",
//     "mlst",
//     "plasmid", 
//     "abritamr",
//     "assembly",
//     "typing",
//     "pangenome"
//     "preview"
//    "tb"
include { READ_ANALYSIS } from './workflows/seq_assessment'
include { RUN_ASSEMBLE } from './workflows/assemble'
workflow {
    // Sequence assessment for reads and assemblies
    // then depending on the workflow run the appropriate analysis
    println "The samples are : $samples"
    reads_pe = Channel.fromFilePairs(["${params.outdir}/*/R{1,2}.fastq.gz"])
                    .filter { sample, files -> samples.containsKey(files[0].getParent().getName())}
                    .map { sample, files -> tuple([id: files[0].getParent().getName(),modules: samples[files[0].getParent().getName()],asm:false, pe_reads:true], files)}
    asm = Channel.fromPath( ["${params.outdir}/*/contigs.fa"])       
                .filter { files -> samples.containsKey(files.getParent().getName())} 
                .map {  files -> tuple([id: files.getParent().getName(), modules:samples[files.getParent().getName()], asm:true, pe_reads:false], files)}
    
    




    // ASSEMBLY_ANALYSIS ( asm )
    // results = READ_ANALYSIS.out.read_stats.concat( ASSEMBLY_ANALYSIS.out.assembly_stats )

    // need a parameter to select the method for speciation
    // speciation 
    // assembly - only if no assembly is provided
    // snps  
        // is default for tree building
    // ska
    // tree

    // typing
    // pangenome
}
