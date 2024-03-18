#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CNVPYTOR_VIEW } from '../../../../../modules/nf-core/cnvpytor/view/main.nf'

workflow test_cnvpytor_view {

    input = [
        [ id:'test'], // meta map
        [file(params.test_data['homo_sapiens']['illumina']['test_pytor'], checkIfExists: true)]
    ]

    bin_sizes = "10000 100000"

    CNVPYTOR_VIEW ( input, bin_sizes, [] )
}

workflow test_cnvpytor_view_tsvout {

    input = [
        [ id:'test'], // meta map
        [file(params.test_data['homo_sapiens']['illumina']['test_pytor'], checkIfExists: true)]
    ]

    output_suffix = "tsv"

    CNVPYTOR_VIEW ( input, [], output_suffix )
}

workflow test_cnvpytor_view_stub {

    input = [
        [ id:'test'], // meta map
        [file(params.test_data['homo_sapiens']['illumina']['test_pytor'], checkIfExists: true)]
    ]

    bin_sizes = []
    output_suffix = []

    CNVPYTOR_VIEW ( input, bin_sizes, output_suffix )
}
