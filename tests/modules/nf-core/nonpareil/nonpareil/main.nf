#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GUNZIP              } from '../../../../../modules/nf-core/gunzip/main.nf'
include { NONPAREIL_NONPAREIL } from '../../../../../modules/nf-core/nonpareil/nonpareil/main.nf'

workflow test_nonpareil_nonpareil {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['candidatus_portiera_aleyrodidarum']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    GUNZIP (input)

    NONPAREIL_NONPAREIL ( GUNZIP.out.gunzip, 'fastq', 'kmer' )
}

workflow test_nonpareil_nonpareil_nogunzip {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['candidatus_portiera_aleyrodidarum']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    NONPAREIL_NONPAREIL ( input, 'fastq', 'kmer' )
}
