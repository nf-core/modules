#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GFASTATS } from '../../../../modules/nf-core/gfastats/main.nf'

workflow test_gfastats_asm_only {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    GFASTATS (
        input,
        '',      // No genome size
        '',      // No target
        [],      // No agp file
        [],      // No include bed file
        [],      // No exclude bed file
        []       // No manipulation instructions file.
    )
}
