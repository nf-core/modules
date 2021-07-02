#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RSEM_CALCULATEEXPRESSION }   from '../../../../software/rsem/calculateexpression/main.nf'   addParams(options: [:])

workflow test_rsem_calculateexpression {
    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
            ]

    RSEM_CALCULATEEXPRESSION ( input )
}
