#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UCSC_GTFTOGENEPRED } from '../../../../../modules/nf-core/ucsc/gtftogenepred/main.nf'

workflow test_ucsc_gtftogenepred {

    input = [
        [ id: 'test' ], // meta map
        [ file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true ) ]
    ]

    UCSC_GTFTOGENEPRED ( input )
}
