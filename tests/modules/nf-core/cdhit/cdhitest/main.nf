#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CDHIT_CDHITEST } from '../../../../../modules/nf-core/cdhit/cdhitest/main.nf'

workflow test_cdhit_cdhitest {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true)
    ]

    CDHIT_CDHITEST ( input )
}
