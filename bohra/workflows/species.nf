#!/usr/bin/env nextflow
include { RUN_KRAKEN } from './../subworkflows/kraken'
include { EXTRACT_SPECIES } from './../modules/extract_species/main'
include { COMBINE_SPECIES_VALS } from './../modules/combine_species_vals/main'
include { COMBINE_SPECIES_REPORT } from './../modules/combine_species/main'
include { CSVTK_CONCAT;CSVTK_UNIQ } from './../modules/csvtk/main'


workflow RUN_SPECIES_READS {
    take:
        sequences
    main:
        
        if ( params.use_kraken2 ) {
            
            RUN_KRAKEN ( sequences )
            species_raw = RUN_KRAKEN.out.species_raw
            species = RUN_KRAKEN.out.species
            species_obs = RUN_KRAKEN.out.species_obs
            version = RUN_KRAKEN.out.version
        } 
        
    emit:
        species_raw  = species_raw
        species = species
        species_obs = species_obs
        version 

}


workflow RUN_SPECIES_ASM {
    take:
        sequences
    main:
        species = sequences.map { cfg, files -> tuple(cfg, "no_results") }
        species_raw = sequences.map { cfg, files -> tuple(cfg, "no_results") }
        species_obs = sequences.map { cfg, files -> tuple(cfg, "no_results") }
        version = Channel.empty()
        if ( params.use_kraken2 ) {
            RUN_KRAKEN ( sequences )
            species_raw = RUN_KRAKEN.out.species_raw
            species = RUN_KRAKEN.out.species
            species_obs = RUN_KRAKEN.out.species_obs
            version = RUN_KRAKEN.out.version
        
        } 

    emit:
        species_raw = species_raw
        species = species
        species_obs = species_obs
        version = version

}


workflow COMBINE_SPECIES {
    take:
        reads_results
        reads_result_file
        asm_results
        asm_result_file
        
    main:
        
        sp_res = reads_results.map { cfg,species -> tuple(cfg.id, cfg.findAll {it.key != 'input_type'}, species.trim() ) }
        asm_res = asm_results.map { cfg,species -> tuple(cfg.id, cfg.findAll {it.key != 'input_type'}, species.trim() ) }
        
        observed_species =  sp_res.join(asm_res, remainder:true)
                                .map( v -> v.findAll { it != null } )
                                .map( v -> { v.size() < 4 ? v[1..-1] + 'no_results': [v[1],v[2],v[-1]]} )
                                .map { cfg, species_reads, species_asm -> tuple(cfg, species_reads ? species_reads : 'no_results', species_asm ? species_asm: 'no_results') }
        
        COMBINE_SPECIES_VALS ( observed_species )
        species_obs = COMBINE_SPECIES_VALS.out.extracted_species
        
        reads_results_files = reads_result_file.map { cfg,files -> tuple(cfg.id, cfg.findAll {it.key != 'input_type'}, files ) }
        
        asm_results_files = asm_result_file.map { cfg,file -> tuple(cfg.id, cfg.findAll {it.key != 'input_type'}, file ) }
        species = reads_results_files.join( asm_results_files, remainder: true )
                    .map( v -> v.findAll { it != null } )
                                .map( v -> { v.size() < 4 ? v[1..-1] + 'no_results': [v[1],v[2],v[-1]]} )
                                .map { cfg, species_reads, species_asm -> tuple(cfg, species_reads ? species_reads : 'no_results', species_asm ? species_asm: 'no_results') }
        
        COMBINE_SPECIES_REPORT ( species )
        species_report = COMBINE_SPECIES_REPORT.out.species_report
        species_stats = species_report.map { cfg, sp -> sp }.collect()
        species_stats = species_stats.map { files -> tuple("speciation", files) }
        
        CSVTK_CONCAT ( species_stats )


    emit:
        species_obs
        species_summary = CSVTK_CONCAT.out.collated
}