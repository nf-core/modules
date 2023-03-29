#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MINDAGAP_MINDAGAP } from '../../../../../modules/nf-core/mindagap/mindagap/main.nf'

workflow test_mindagap {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['imaging']['tiff']['mouse_heart_wga'], checkIfExists: true)
    ]

    boxsize = 3
    loopnum = 40

    MINDAGAP_MINDAGAP ( input, boxsize, loopnum )
}
