#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MAGECK_TEST } from '../../../../../modules/nf-core/mageck/test/main.nf'

workflow test_mageck_test {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['mus_musculus']['csv']['count_table'], checkIfExists: true)
    ]


    MAGECK_TEST ( input )
}

workflow test_mageck_test_day0_label {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['mus_musculus']['csv']['count_table'], checkIfExists: true)
    ]


    MAGECK_TEST ( input )
}
