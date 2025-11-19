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
def ctrl = [:]
input_file = file(params.isolates) // need to make this an input file 
    reader =   input_file.newReader()
    reader.eachLine { line ->
    if (!line.startsWith("Isolate") & !line.startsWith("#")) {
        
    line = line.split("\t")
        // println "The line is : $line"
        samples[line[0]] = line[0]
        asms[line[0]] = line[3]
        // asmblr[line[0]] = line[3]
        sp[line[0]] = line[4]
        ctrl[line[0]] = line[5] 
     }
    }

println "Will run : $params.modules"
// println samples
include { READ_ANALYSIS;ASSEMBLY_ANALYSIS } from './workflows/seq_assessment'
include { RUN_ASSEMBLE } from './workflows/assemble'
include { RUN_SPECIES_READS; RUN_SPECIES_ASM; COMBINE_SPECIES } from './workflows/species'
include { RUN_TYPING } from './workflows/typing'
include { RUN_TBTAMR } from './workflows/tbtamr'
include { RELATIONSHIPS } from './workflows/relationships'
include { RUN_PANAROO } from './workflows/pangenome'
include { RUN_COMPILE } from './workflows/compile'

workflow {
    // Sequence assessment for reads and assemblies
    // then depending on the workflow run the appropriate analysis
    reads_pe = Channel.fromFilePairs(["${params.outdir}/*/R{1,2}.fastq.gz"])
                    .filter { sample, files -> samples.containsKey(files[0].getParent().getName())}
                    .map { sample, files -> tuple([id: files[0].getParent().getName(),input_type:'pe_reads', asm :asms[files[0].getParent().getName()],  species:sp[files[0].getParent().getName()], control:ctrl[files[0].getParent().getName()]], files)}
    // reads_ont = Channel.fromPath( ["${params.outdir}/*/read_ont.fastq.gz","${params.outdir}/*/read_ont.fastq"])       
    //             .filter { files -> samples.containsKey(files.getParent().getName())} 
    //             .map {  files -> tuple([id: files.getParent().getName(), modules:samples[files.getParent().getName()], input_type:'ont_reads', asm :asms[files.getParent().getName()],assembler:asmblr[files.getParent().getName()],species:sp[files.getParent().getName()]], files)}
    
    asm = Channel.fromPath( ["${params.outdir}/*/contigs.fa"])       
                .filter { files -> samples.containsKey(files.getParent().getName())} 
                .map {  files -> tuple([id: files.getParent().getName(), input_type:'asm', asm :files, species:sp[files.getParent().getName()], control: ctrl[files.getParent().getName()]], files)}
    
    
    results = Channel.empty()

    // println reads_pe.view()
    // println asm.view()
    // println reads_pe.view()
    // seq assessment is always done on every input
    READ_ANALYSIS ( reads_pe )
    read_stats = READ_ANALYSIS.out.read_stats
    reads_pe = READ_ANALYSIS.out.reads_pe
    results = results.concat( read_stats)
    versions = READ_ANALYSIS.out.version_seqkit_reads
    versions = versions.concat( READ_ANALYSIS.out.version_kmc )
    versions = versions.concat( READ_ANALYSIS.out.version_bohra )
    // reads_pe = READ_ANALYSIS.out.reads_pe
    // if there is assembly in the modules list then generate an assembly and run assembly analysis
    
    if (params.modules.contains("assemble") ){
        // assembly is only done if the input is reads
        RUN_ASSEMBLE ( reads_pe.filter { cfg, files -> cfg.control != 'control' } )
        // RUN_ASSEMBLE ( reads_ont )
        asm = RUN_ASSEMBLE.out.contigs
        versions = versions.concat( RUN_ASSEMBLE.out.versions )
        results = results.concat( RUN_ASSEMBLE.out.insertiqr )
        
    } 
    ASSEMBLY_ANALYSIS ( asm )
    assembly_stats = ASSEMBLY_ANALYSIS.out.assembly_stats
    versions = versions.concat( ASSEMBLY_ANALYSIS.out.version_prokka, ASSEMBLY_ANALYSIS.out.version_seqkit_asm )
    // update the results with the assembly stats
    results = results.concat( assembly_stats )
    
    if (params.modules.contains("species") ){
        
        RUN_SPECIES_READS ( reads_pe )
        RUN_SPECIES_ASM ( asm )
        versions = versions.concat( RUN_SPECIES_READS.out.version, RUN_SPECIES_ASM.out.version )
        reads_species_obs = RUN_SPECIES_READS.out.species_obs
        asm_species_obs = RUN_SPECIES_ASM.out.species_obs
        sp_res = reads_species_obs.map { cfg,species -> tuple(cfg.id, cfg.findAll {it.key != 'input_type'}, species.trim() ) }
        asm_res = asm_species_obs.map { cfg,species -> tuple(cfg.id, cfg.findAll {it.key != 'input_type'}, species.trim() ) }
        COMBINE_SPECIES ( 
            RUN_SPECIES_READS.out.species_obs, 
            RUN_SPECIES_READS.out.species,
            RUN_SPECIES_ASM.out.species_obs, 
            RUN_SPECIES_ASM.out.species
            )
        
       
        species_tmp = COMBINE_SPECIES.out.species_obs
                                            .map { cfg, species_obs -> tuple(cfg.id, cfg, species_obs.trim() ) }
        // println species_tmp.view()
        reads = reads_pe.map { cfg, files -> tuple(cfg.id, cfg, files) }
        reads_pe = reads.join( species_tmp )
                                .map { id, cfg_reads, files, cfg_spieces,  species_obs -> tuple(cfg_reads + [species:species_obs.trim()] , files) }
                         
        asm_tmp = asm.map { cfg, files -> tuple(cfg.id, cfg , files) }
        asm = asm_tmp.join( species_tmp )
                                .map { id, cfg_asm, files, cfg_spieces, species_obs -> tuple(cfg_asm + [species:species_obs.trim()] , files) }
        
        // // generate summay file for species
        species_report = COMBINE_SPECIES.out.species_summary
        results = results.concat( species_report )
        
    }
    
    if (params.modules.contains("mtb")){
        RUN_TBTAMR ( reads_pe.filter { cfg, files -> cfg.control != 'control' } )
        results = results.concat( RUN_TBTAMR.out.results )
        versions = versions.concat( RUN_TBTAMR.out.version )
    }

    if (params.modules.contains("typing") ){
        // assembly is only done if the input is reads
        // find any tb as plasmid and mlst and abritamr no good - use tbtamr
        asm_typing = asm.filter { cfg, asm -> cfg.species != 'Mycobacterium tuberculosis' }.filter { cfg, files -> cfg.control != 'control' }
        reads_nottb = reads_pe.filter { cfg, reads -> cfg.species != 'Mycobacterium tuberculosis' }.filter { cfg, files -> cfg.control != 'control' }
        reads_tb = reads_pe.filter { cfg, reads -> cfg.species == 'Mycobacterium tuberculosis' }.filter { cfg, files -> cfg.control != 'control' }
        RUN_TYPING ( asm_typing, reads_nottb )
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
        versions = versions.concat( RUN_TYPING.out.versions )
        RUN_TBTAMR ( reads_tb)
        results = results.concat( RUN_TBTAMR.out.results )
        versions = versions.concat( RUN_TBTAMR.out.version )
        
    }
    
    if (params.modules.contains("snippy") || (params.modules.contains("ska")) || (params.modules.contains("mash"))){
        
        if (params.modules.contains("snippy") ){
            sequences = reads_pe.filter { cfg, files -> cfg.control != 'control' }
            
            // RUN_SNIPPY ( reads_pe, Channel.from(file(params.reference)) )
        } else if (params.modules.contains("ska") || params.modules.contains("mash") ){
            reads = reads_pe.filter { cfg,files -> cfg.control != 'control' }.map { cfg, files -> tuple(cfg.id, cfg, files) }
            // 
            asm_tmp = asm.filter { cfg,files -> cfg.control != 'control' }.map { cfg, files -> tuple(cfg.id, cfg , files) }
            // println reads.view()
            sequences = reads.join(asm_tmp, remainder:true).map( v -> { v.size() == 4 ? v[1] ? [v[1],v[2]] : [v[2],v[3]]  : [v[1],v[2]]} )
           
        } 
        reference = Channel.fromPath(params.reference).ifEmpty('no_file')
        // println sequences.view()
        RELATIONSHIPS ( sequences, reference )
      
        results = results.concat( RELATIONSHIPS.out.dists )
        results = results.concat( RELATIONSHIPS.out.core_vcf )
        results = results.concat( RELATIONSHIPS.out.clusters )
        results = results.concat( RELATIONSHIPS.out.stats )
        results = results.concat( RELATIONSHIPS.out.tree)
        versions = versions.concat( RELATIONSHIPS.out.version )
        versions = versions.concat( RELATIONSHIPS.out.tree_version )

    }

    if (params.modules.contains("pangenome")){

        if( params.pangenome_groups == "clusters") {
            groups = RELATIONSHIPS.out.clusters
        }
        gff = ASSEMBLY_ANALYSIS.out.gff
        RUN_PANAROO ( gff,groups )

        
        results = results.concat( RUN_PANAROO.out.pangenome_rtab )
        results = results.concat( RUN_PANAROO.out.classification )
        results = results.concat( RUN_PANAROO.out.groups )
        versions = versions.concat( RUN_PANAROO.out.version )
    }

    // println results.view()
    RUN_COMPILE ( results, versions )
}
