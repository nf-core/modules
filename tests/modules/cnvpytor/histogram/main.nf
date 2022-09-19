#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CNVPYTOR_HISTOGRAM } from '../../../../modules/cnvpytor/histogram/main.nf'

workflow test_cnvpytor_histogram {

    input = [
        [ id:'test'], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_pytor'], checkIfExists: true),
    ]

    CNVPYTOR_HISTOGRAM ( input )
}
