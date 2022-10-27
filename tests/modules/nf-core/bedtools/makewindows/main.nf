#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { BEDTOOLS_MAKEWINDOWS } from "$moduleDir/modules/nf-core/bedtools/makewindows/main.nf"

workflow test_bedtools_makewindows {

    input = [
        [ id:'test'],
        file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
    ]

    BEDTOOLS_MAKEWINDOWS ( input, true )
}
