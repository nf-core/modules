#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MARKDUP_PARALLEL } from '../../../../subworkflows/nf-core/markdup_parallel/main' addParams([:])

workflow test_markdup_parallel {
    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
            ]

    MARKDUP_PARALLEL ( input)
}
