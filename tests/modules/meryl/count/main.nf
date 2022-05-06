#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MERYL_COUNT } from '../../../../modules/meryl/count/main.nf'

workflow test_meryl_count_single_end {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    MERYL_COUNT ( input )
}
