#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PHISPY } from '../../../../modules/nf-core/phispy/main.nf'

workflow test_phispy {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['bacteroides_fragilis']['genome']['genome_gbff_gz'], checkIfExists: true)
    ]

    PHISPY ( input )
}
