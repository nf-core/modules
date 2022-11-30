#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ENTREZDIRECT_ESUMMARY } from '../../../../../modules/nf-core/entrezdirect/esummary/main.nf'

//
// Test with SRA database
//
workflow test_entrezdirect_esummary_sra {

    input = [
        [ id:'test_sra' ], // meta map
        uid = '5135484',
        []
    ]
    database = 'sra'

    ENTREZDIRECT_ESUMMARY ( input, database )
}

//
// Test with Genome database
//
workflow test_entrezdirect_esummary_genome {

    input = [
        [ id:'test_genome' ], // meta map
        uid = '768',
        []
    ]
    database = 'genome'

    ENTREZDIRECT_ESUMMARY ( input, database )
}

//
// Test with Assembly database
//
workflow test_entrezdirect_esummary_assembly {

    input = [
        [ id:'test_assembly' ], // meta map
        uid = '191021',
        []
    ]
    database = 'assembly'

    ENTREZDIRECT_ESUMMARY ( input, database )
}
