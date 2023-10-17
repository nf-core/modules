#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { WISECONDORX_NEWREF  } from '../../../../../modules/nf-core/wisecondorx/newref/main.nf'
include { WISECONDORX_CONVERT } from '../../../../../modules/nf-core/wisecondorx/convert/main.nf'

workflow test_wisecondorx_newref {
    
    input = Channel.of([
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ],
    [
        [ id: 'test2' ],
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true)
    ])

    WISECONDORX_CONVERT(
        input,
        [[],[]],
        [[],[]]
    )

    WISECONDORX_CONVERT.out.npz
        .map { meta, npz ->
            new_meta = meta + [id:"combined"]
            [new_meta, npz]
        }
        .groupTuple()
        .set { newref_input }

    WISECONDORX_NEWREF ( newref_input )
}
