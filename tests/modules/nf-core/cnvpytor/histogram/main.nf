#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { CNVPYTOR_HISTOGRAM } from "$moduleDir/modules/nf-core/cnvpytor/histogram/main.nf"

workflow test_cnvpytor_histogram {

    input = [
        [ id:'test'], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_pytor'], checkIfExists: true),
    ]

    CNVPYTOR_HISTOGRAM ( input )
}
