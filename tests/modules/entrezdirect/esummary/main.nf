#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ENTREZDIRECT_ESUMMARY } from '../../../../modules/entrezdirect/esummary/main.nf'

//
// Test with SRA database
//
workflow test_entrezdirect_esummary_sra {

    input_file = file('NO_FILE')

    input = [
        [ id:'test_sra' ], // meta map
        uid = '5135484',
        uids_file = input_file
    ]

    database = 'sra'

    ENTREZDIRECT_ESUMMARY ( input, database )
}

//
// Test with Genome database
//
workflow test_entrezdirect_esummary_genome {

    input_file = file('NO_FILE')

    input = [
        [ id:'test_genome' ], // meta map
        uid = '768',
        uids_file = input_file
    ]

    database = 'genome'

    ENTREZDIRECT_ESUMMARY ( input, database )
}

//
// Test with Assembly database
//
workflow test_entrezdirect_esummary_assembly {

    input_file = file('NO_FILE')

    input = [
        [ id:'test_assembly' ], // meta map
        uid = '191021',
        uids_file = input_file
    ]

    database = 'assembly'

    ENTREZDIRECT_ESUMMARY ( input, database )
}
