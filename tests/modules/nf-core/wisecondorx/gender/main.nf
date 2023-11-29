#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { WISECONDORX_GENDER  } from '../../../../../modules/nf-core/wisecondorx/gender/main.nf'
include { WISECONDORX_CONVERT } from '../../../../../modules/nf-core/wisecondorx/convert/main.nf'

workflow test_wisecondorx_gender {
    
    input = Channel.of([
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam_bai'], checkIfExists: true)
    ])

    WISECONDORX_CONVERT(
        input,
        [[], []],
        [[], []]
    )

    reference = [
        [ id: "reference" ],
        file("reference.npz")
    ]

    WISECONDORX_GENDER(
        WISECONDORX_CONVERT.out.npz,
        reference
    )
}
