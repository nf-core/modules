#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NONPAREIL_NONPAREIL } from '../../../../../modules/nf-core/nonpareil/nonpareil/main.nf'
include { NONPAREIL_CURVE } from '../../../../../modules/nf-core/nonpareil/curve/main.nf'

workflow test_nonpareil_curve {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['candidatus_portiera_aleyrodidarum']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    NONPAREIL_NONPAREIL ( input, 'fastq', 'kmer' )
    NONPAREIL_CURVE (NONPAREIL_NONPAREIL.out.npo)

}
