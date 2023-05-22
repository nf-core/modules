#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PURECN_COVERAGE } from '../../../../../modules/nf-core/purecn/coverage/main.nf'

workflow test_purecn_coverage {

    meta = [ id:'test' ]
    input_bam = [ meta,
                  file(params.test_data['homo_sapiens']['illumina']['purecn_ex1_bam'], checkIfExists: true),
                  file(params.test_data['homo_sapiens']['illumina']['purecn_ex1_bai'], checkIfExists: true)
                  ]
    input_intervals = file(params.test_data['homo_sapiens']['illumina']['purecn_ex1_interval'], checkIfExists: true)


    PURECN_COVERAGE ( input_bam, input_intervals )
}
