#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MEGAN_RMA2INFO } from '../../../../../modules/nf-core/megan/rma2info/main.nf'

workflow test_megan_rma2info {

    input = [
        [ id:'test', single_end:false ], // meta map
        file('https://github.com/nf-core/test-datasets/raw/a7e61654553887475a2f7178108587ecd9b54608/data/delete_me/malt/test.rma6', checkIfExists: true)
    ]
    megan_summary = true

    MEGAN_RMA2INFO ( input, megan_summary )
}
