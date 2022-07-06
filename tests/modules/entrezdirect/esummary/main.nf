#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ENTREZDIRECT_ESUMMARY } from '../../../../modules/entrezdirect/esummary/main.nf'

//
// Test with SRA database
//
workflow test_entrezdirect_esummary_sra {
    input = [
        [ id:'test_sra' ], // meta map
        database = 'sra',
        uid = '5135484'
    ]

    ENTREZDIRECT_ESUMMARY ( input )
}

//
// Test with Genome database
//
workflow test_entrezdirect_esummary_genome {
    input = [
        [ id:'test_genome' ], // meta map
        database = 'genome',
        uid = '768'
    ]

    ENTREZDIRECT_ESUMMARY ( input )
}

//
// Test with Assembly database
//
workflow test_entrezdirect_esummary_assembly {
    input = [
        [ id:'test_assembly' ], // meta map
        database = 'assembly',
        uid = '191021'
    ]

    ENTREZDIRECT_ESUMMARY ( input )
}
