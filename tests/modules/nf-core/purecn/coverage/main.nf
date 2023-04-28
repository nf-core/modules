#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PURECN_COVERAGE } from '../../../../../modules/nf-core/purecn/coverage/main.nf'

workflow test_purecn_coverage {

    meta = [ id:'test' ]
    input_bam = [ meta,
                  file("https://raw.githubusercontent.com/lima1/PureCN/master/inst/extdata/ex1.bam", checkIfExists: true),
                  file("https://raw.githubusercontent.com/lima1/PureCN/master/inst/extdata/ex1.bam.bai", checkIfExists: true)
                  ]
    input_intervals = ["https://raw.githubusercontent.com/lima1/PureCN/master/inst/extdata/ex1_intervals.txt"]


    PURECN_COVERAGE ( input_bam, input_intervals )
}
