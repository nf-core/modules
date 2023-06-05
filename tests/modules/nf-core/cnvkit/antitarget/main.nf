#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CNVKIT_ANTITARGET } from '../../../../../modules/nf-core/cnvkit/antitarget/main.nf'

workflow test_cnvkit_antitarget {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true)
    ]

    CNVKIT_ANTITARGET ( input )
}

