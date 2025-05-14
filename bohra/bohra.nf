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
def asms = [:]
def asmblr = [:]
def sp = [:]
input_file = file(params.isolates) // need to make this an input file 
    reader =   input_file.newReader()
    reader.eachLine { line ->
    line = line.split("\t")
        // println "The line is : $line"
        samples[line[0]] = line[1]
        asms[line[0]] = line[2]
        asmblr[line[0]] = line[3]
        sp[line[0]] = line[4]
    }
// get reads
println "Will run : $params.modules"
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
include { READ_ANALYSIS;ASSEMBLY_ANALYSIS } from './workflows/seq_assessment'
include { RUN_ASSEMBLE } from './workflows/assemble'
include { RUN_SPECIES_READS; RUN_SPECIES_ASM; COMBINE_SPECIES } from './workflows/species'
include { RUN_TYPING } from './workflows/typing'
include { RELATIONSHIPS } from './workflows/relationships'
include { PREVIEW_NEWICK } from './workflows/preview'

workflow {
    // Sequence assessment for reads and assemblies
    // then depending on the workflow run the appropriate analysis
    // println "The samples are : $samples"
    reads_pe = Channel.fromFilePairs(["${params.outdir}/*/R{1,2}.fastq.gz"])
                    .filter { sample, files -> samples.containsKey(files[0].getParent().getName())}
                    .map { sample, files -> tuple([id: files[0].getParent().getName(),modules: samples[files[0].getParent().getName()],input_type:'pe_reads', asm :asms[files[0].getParent().getName()], assembler:asmblr[files[0].getParent().getName()], species:sp[files[0].getParent().getName()] ], files)}
    reads_ont = Channel.fromPath( ["${params.outdir}/*/read_ont.fastq.gz","${params.outdir}/*/read_ont.fastq"])       
                .filter { files -> samples.containsKey(files.getParent().getName())} 
                .map {  files -> tuple([id: files.getParent().getName(), modules:samples[files.getParent().getName()], input_type:'ont_reads', asm :asms[files.getParent().getName()],assembler:asmblr[files.getParent().getName()],species:sp[files.getParent().getName()]], files)}
    
    asm = Channel.fromPath( ["${params.outdir}/*/contigs.fa"])       
                .filter { files -> samples.containsKey(files.getParent().getName())} 
                .map {  files -> tuple([id: files.getParent().getName(), modules:samples[files.getParent().getName()], input_type:'asm', asm :files,assembler:asmblr[files.getParent().getName()],species:sp[files.getParent().getName()]], files)}
    
    
    results = Channel.empty()
    // seq assessment is always done on every input
    READ_ANALYSIS ( reads_pe )
    read_stats = READ_ANALYSIS.out.read_stats
    results = results.concat( read_stats)
    // if there is assembly in the modules list then generate an assembly and run assembly analysis
    if (params.modules.contains("assemble") ){
        // assembly is only done if the input is reads
        RUN_ASSEMBLE ( reads_pe )
        // RUN_ASSEMBLE ( reads_ont )
        asm = RUN_ASSEMBLE.out.contigs
        
    } 
    ASSEMBLY_ANALYSIS ( asm )
    assembly_stats = ASSEMBLY_ANALYSIS.out.assembly_stats
    // update the results with the assembly stats
    results = results.concat( assembly_stats )

    if (params.modules.contains("species") ){
        
        RUN_SPECIES_READS ( reads_pe )
        RUN_SPECIES_ASM ( asm )
        
        reads_species_obs = RUN_SPECIES_READS.out.species_obs
        asm_species_obs = RUN_SPECIES_ASM.out.species_obs
        // println asm_species_obs.view()
        // println reads_species_obs.view()
        sp_res = reads_species_obs.map { cfg,species -> tuple(cfg.id, cfg.findAll {it.key != 'input_type'}, species.trim() ) }
        asm_res = asm_species_obs.map { cfg,species -> tuple(cfg.id, cfg.findAll {it.key != 'input_type'}, species.trim() ) }
        // println sp_res.join(asm_res, remainder:true).map( v -> { v.size() == 4 ? v[1..3] : [v[1],v[2],v[4]]} ).view()
        COMBINE_SPECIES ( 
            RUN_SPECIES_READS.out.species_obs, 
            RUN_SPECIES_READS.out.species,
            RUN_SPECIES_ASM.out.species_obs, 
            RUN_SPECIES_ASM.out.species
            )
        
       
        species_tmp = COMBINE_SPECIES.out.species_obs
                                            .map { cfg, species_obs -> tuple(cfg.id, cfg, species_obs.trim() ) }
        reads = reads_pe.map { cfg, files -> tuple(cfg.id, cfg, files) }
        reads_pe = reads.join( species_tmp )
                                .map { id, cfg_reads, files, cfg_spieces,  species_obs -> tuple(cfg_reads + [species:species_obs.trim()] , files) }
                            
        asm_tmp = asm.map { cfg, files -> tuple(cfg.id, cfg , files) }
        asm = asm_tmp.join( species_tmp )
                                .map { id, cfg_asm, files, cfg_spieces, species_obs -> tuple(cfg_asm + [species:species_obs.trim()] , files) }
        
        // generate summay file for species
        species_report = COMBINE_SPECIES.out.species_summary
        results = results.concat( species_report )
        
    }
    
    if (params.modules.contains("typing") ){
        // assembly is only done if the input is reads
        RUN_TYPING ( asm, reads_pe )
        resistome = RUN_TYPING.out.resistome
        virulome = RUN_TYPING.out.virulome
        plasmid = RUN_TYPING.out.plasmid
        inferred = RUN_TYPING.out.inferred.ifEmpty { "no_results" }
        reportable = RUN_TYPING.out.reportable
        serotypes = RUN_TYPING.out.serotypes
        mlst = RUN_TYPING.out.mlst
        results = results.concat( resistome )
        results = results.concat( virulome )
        results = results.concat( plasmid )
        results = results.concat( inferred )
        results = results.concat( reportable )
        results = results.concat( serotypes )
        results = results.concat( mlst )
    }
    if (params.modules.contains("mash")) {
        reads = reads_pe.map { cfg, files -> tuple(cfg.id, cfg, files) }
        asm_tmp = asm.map { cfg, files -> tuple(cfg.id, cfg , files) }
        sequences = reads.join(asm_tmp, remainder:true).map( v -> { v.size() == 4 ? v[1] ? [v[1],v[2]] : [v[2],v[3]]  : [v[1],v[2]]} )
        PREVIEW_NEWICK ( sequences )
        results = results.concat( PREVIEW_NEWICK.out.nwk )
    }
    if (params.modules.contains("snippy") || (params.modules.contains("ska"))){
        
        if (params.modules.contains("snippy") ){
            sequences = reads_pe
            // RUN_SNIPPY ( reads_pe, Channel.from(file(params.reference)) )
        } else if (params.modules.contains("ska") ){
            reads = reads_pe.map { cfg, files -> tuple(cfg.id, cfg, files) }
            asm_tmp = asm.map { cfg, files -> tuple(cfg.id, cfg , files) }
            sequences = reads.join(asm_tmp, remainder:true).map( v -> { v.size() == 4 ? v[1] ? [v[1],v[2]] : [v[2],v[3]]  : [v[1],v[2]]} )
           
        }
        RELATIONSHIPS ( sequences, Channel.fromPath(params.reference) )
        // reads = reads_pe.map { cfg, files -> tuple(cfg.id, cfg, files) }
        // asm_tmp = asm.map { cfg, files -> tuple(cfg.id, cfg , files) }

        // sp_res = reads_results.map { cfg,species -> tuple(cfg.id, cfg.findAll {it.key != 'input_type'}, species.trim() ) }
        // asm_res = asm_results.map { cfg,species -> tuple(cfg.id, cfg.findAll {it.key != 'input_type'}, species.trim() ) }
        
        // observed_species =  sp_res.join(asm_res, remainder:true).map( v -> { v.size() == 4 ? v[1..3] : [v[1],v[2],v[4]]} )
        //                                                          .map { cfg, species_reads, species_asm -> tuple(cfg, species_reads ? species_reads : 'no_results', species_asm ? species_asm: 'no_results') }
        

    }

    // snps  
        // is default for tree building
    // ska
    // tree

    // typing
    // pangenome
}
