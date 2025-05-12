#!/usr/bin/env nextflow

workflow ABRITAMR {
    take:
        input
    main:
        ABRITAMR ( input )
    emit:
        resistomes = ABRITAMR.out.resistomes
        virulomes = ABRITAMR.out.virulomes
}