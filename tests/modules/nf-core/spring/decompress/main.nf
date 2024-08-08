#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SPRING_COMPRESS } from '../../../../../modules/nf-core/spring/compress/main.nf'
include { SPRING_DECOMPRESS } from '../../../../../modules/nf-core/spring/decompress/main.nf'

workflow test_spring_decompress_single_end {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
        []
    ]

    SPRING_COMPRESS ( input )
    SPRING_DECOMPRESS ( SPRING_COMPRESS.out.spring, true )
}

workflow test_spring_decompress_paired_end {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
    ]

    SPRING_COMPRESS ( input )
    SPRING_DECOMPRESS ( SPRING_COMPRESS.out.spring, false )
}
