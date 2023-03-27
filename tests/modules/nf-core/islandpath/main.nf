#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ISLANDPATH } from '../../../../modules/nf-core/islandpath/main.nf'

workflow test_islandpath {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['bacteroides_fragilis']['genome']['genome_gbff_gz'], checkIfExists: true)
    ]

    ISLANDPATH ( input )
}
