#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { DEEPTOOLS_COMPUTEMATRIX } from "$moduleDir/modules/nf-core/deeptools/computematrix/main.nf"

workflow test_deeptools_computematrix {

    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_bigwig'], checkIfExists: true)
            ]

    bed = file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)

    DEEPTOOLS_COMPUTEMATRIX ( input, bed )
}
