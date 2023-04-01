#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PURECN_COVERAGE } from '../../../../../modules/nf-core/purecn/coverage/main.nf'

workflow test_purecn_coverage {

    meta = [ id:'test' ]
    input_bam = [ meta, file(params.test_data['homo_sapienza']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ]
    input_intervals = []//intervals file from purecn/intervalfile need to be uploaded in the test data repo


    PURECN_COVERAGE ( input_bam, input_intervals )
}
