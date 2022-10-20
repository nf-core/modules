#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MERYL_COUNT     } from '../../../../../modules/nf-core/meryl/count/main.nf'
include { MERYL_HISTOGRAM } from '../../../../../modules/nf-core/meryl/histogram/main.nf'

workflow test_meryl_histogram {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    MERYL_COUNT ( input )
    MERYL_HISTOGRAM ( MERYL_COUNT.out.meryl_db )
}
