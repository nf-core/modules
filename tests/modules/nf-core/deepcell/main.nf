#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DEEPCELL } from '../../../../modules/nf-core/deepcell/main.nf'

workflow test_deepcell {

    img = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['image'], checkIfExists: true)
    ]

    DEEPCELL ( img )
}
