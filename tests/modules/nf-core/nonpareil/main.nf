#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GUNZIP } from '../../../../modules/nf-core/gunzip/main.nf'
include { NONPAREIL } from '../../../../modules/nf-core/nonpareil/main.nf'

workflow test_nonpareil {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['candidatus_portiera_aleyrodidarum']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    GUNZIP (input)

    NONPAREIL ( GUNZIP.out.gunzip, 'fastq', 'kmer' )
}
