#!/usr/bin/env nextflow
include { ABRITAMR; ABRITAMR_GENERAL;ABRITAMR_INFER;COMBINE_AMR } from './../modules/abritamr/main'
include { JSON_COMBINE } from './../modules/json_combine/main'
include { CSVTK_CONCAT; CSVTK_UNIQ } from './../modules/csvtk/main'
include { CONCAT_FILES } from './../modules/utils/main'


workflow RUN_ABRITAMR {
    
    take:
        asm
        plamsid
        
    main:
        
        ABRITAMR ( asm )
        abritamr_matches = ABRITAMR.out.abritamr_matches
                                            .map { cfg, file -> tuple(cfg.id, cfg, file) }
        abritamr_partials = ABRITAMR.out.abritamr_partials
                                            .map { cfg, file -> tuple(cfg.id, cfg, file) }
        abritamr_summary = abritamr_matches.join( abritamr_partials)
                                            .map {id, cfg_matches, file_mathces, cfg_partials, file_partials-> tuple(cfg_matches, file_mathces, file_partials) }  
        for_inferrence = abritamr_matches
                                    .map {id,  cfg, file -> tuple(cfg, file) }
                                    .filter {cfg, file -> cfg.species =~ 'Salmonella' }
                                            
        ABRITAMR_GENERAL ( abritamr_summary )
        ABRITAMR_INFER ( for_inferrence )
        
        
        
        plasmid = plamsid.map { cfg, file -> tuple(cfg.id, cfg, file) }
        amrout = ABRITAMR.out.amrfinder_out.map { cfg, file -> tuple(cfg.id, cfg, file) }

        for_collation = abritamr_matches.join ( abritamr_partials )
                                            .join ( amrout )
                                            .join ( plasmid )
                                            .map {id, cfg_matches, file_mathces, cfg_partials, file_partials, cfg_plasmid, file_plasmid, cfg_amrout, file_amrout -> tuple(cfg_matches, file_mathces, file_partials, file_plasmid, file_amrout) }
                                            
        COMBINE_AMR ( for_collation )
        versions = ABRITAMR.out.version.map { cfg, file -> file }.collect()
                                        .map { files -> tuple("version_abritamr", files) }
        CSVTK_UNIQ ( versions )
        // println ABRITAMR.out.amrfinder_out.view()
    emit:
        
        abritamr_virulence = ABRITAMR.out.abritamr_virulence
        abritamr_reportable= ABRITAMR_GENERAL.out.reportable
        abritamr_infer = ABRITAMR_INFER.out.inferred
        plasmid = COMBINE_AMR.out.plasmid
        resistome = COMBINE_AMR.out.resistome
        version = CSVTK_UNIQ.out.collated
}

workflow CONCAT_RESISTOMES {
    take:
        resistomes
    main:
        CONCAT_FILES ( resistomes )
    emit:
        collated_resistomes = CONCAT_FILES.out.collated
}

workflow CONCAT_VIRULENCE {
    take:
        virulence
    main:
        CONCAT_FILES ( virulence )
    emit:
        collated_virulome = CONCAT_FILES.out.collated
}

workflow CONCAT_PLASMID {
    take:
        plasmids
    main:
        JSON_COMBINE ( plasmids )
    emit:
        collated_plasmid = JSON_COMBINE.out.collated
}

workflow CONCAT_REPORTABLE {
    take:
        reportables
    main:
        CONCAT_FILES ( reportables )
    emit:
        collated_reportable = CONCAT_FILES.out.collated
}

workflow CONCAT_INFERRED {
    take:
        inferred
    main:
        CONCAT_FILES ( inferred )
    emit:
        collated_inferred = CONCAT_FILES.out.collated
}