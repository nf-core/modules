#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_SPLITCRAM } from '../../../../../modules/nf-core/gatk4/splitcram/main.nf'

workflow test_gatk4_splitcram {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true)
    ]

    GATK4_SPLITCRAM ( input )
}
