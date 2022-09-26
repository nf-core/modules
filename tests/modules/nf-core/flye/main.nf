#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FLYE } from '../../../../modules/nf-core/flye/main.nf'

workflow test_flye_pacbio_raw {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['pacbio']['hifi'], checkIfExists: true)
    ]
    mode = "--pacbio-raw"

    FLYE ( input, mode )
}

workflow test_flye_pacbio_corr {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['pacbio']['hifi'], checkIfExists: true)
    ]
    mode = "--pacbio-corr"

    FLYE ( input, mode )
}

workflow test_flye_pacbio_hifi {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['pacbio']['hifi'], checkIfExists: true)
    ]
    mode = "--pacbio-hifi"

    FLYE ( input, mode )
}

workflow test_flye_nano_raw {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['pacbio']['hifi'], checkIfExists: true)
    ]
    mode = "--nano-raw"

    FLYE ( input, mode )
}

workflow test_flye_nano_corr {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['pacbio']['hifi'], checkIfExists: true)
    ]
    mode = "--nano-corr"

    FLYE ( input, mode )
}

workflow test_flye_nano_hq {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['pacbio']['hifi'], checkIfExists: true)
    ]
    mode = "--nano-hq"

    FLYE ( input, mode )
}
