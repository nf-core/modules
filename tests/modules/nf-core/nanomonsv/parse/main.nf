#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NANOMONSV_PARSE } from '../../../../../modules/nf-core/nanomonsv/parse/main.nf'

workflow test_nanomonsv_parse {

    def data_path = 'https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/nanomonsv'
    input = channel.of(
    [
        [ id:'control' ], // meta map
        file("${data_path}/input_bam/test_ctrl.bam", checkIfExists: true),
        file("${data_path}/input_bam/test_ctrl.bam.bai", checkIfExists: true)
    ],
    [
        [ id:'tumor' ], // meta map
        file("${data_path}/input_bam/test_tumor.bam", checkIfExists: true),
        file("${data_path}/input_bam/test_tumor.bam.bai", checkIfExists: true)
    ]
    )

    NANOMONSV_PARSE ( input )
}
