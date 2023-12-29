#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ESTSFS } from '../../../../modules/nf-core/estsfs/main.nf'

workflow test_estsfs {
    
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['generic']['estsfs']['config_file'], checkIfExists: true), file(params.test_data['generic']['estsfs']['data_file'], checkIfExists: true), file(params.test_data['generic']['estsfs']['seed_file'], checkIfExists: true)
    ]

    ESTSFS ( input )
}
