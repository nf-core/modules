#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { WISECONDORX_PREDICT } from '../../../../../modules/nf-core/wisecondorx/predict/main.nf'
include { WISECONDORX_CONVERT } from '../../../../../modules/nf-core/wisecondorx/convert/main.nf'

workflow test_wisecondorx_predict {
    
    sample = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam_bai'], checkIfExists: true)
    ]

    WISECONDORX_CONVERT (
        sample,
        [[],[]],
        [[],[]]
    )

    reference = [
        [ id:'reference' ],
        file("reference.npz")
    ]

    WISECONDORX_PREDICT (
        WISECONDORX_CONVERT.out.npz,
        reference,
        [[],[]]
    )
}
