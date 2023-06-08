#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQKIT_GREP } from '../../../../../modules/nf-core/seqkit/grep/main.nf'

workflow test_seqkit_grep {

    pattern = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_header'], checkIfExists: true)
    ]

    reference = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    ]

    SEQKIT_GREP ( pattern, reference)
}
