#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CNVPYTOR_HISTOGRAM } from '../../../../modules/cnvpytor/histogram/main.nf'

workflow test_cnvpytor_histogram {

    input = [
        [ id:'test'], // meta map
        file(params.pytor_file, checkIfExists: true)
    ]

    CNVPYTOR_HISTOGRAM ( input )
}
