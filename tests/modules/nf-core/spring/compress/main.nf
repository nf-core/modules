#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SPRING_COMPRESS } from '../../../../../modules/nf-core/spring/compress/main.nf'

workflow test_spring_compress_single_end {

    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
        []
    ]

    SPRING_COMPRESS ( input )
}

workflow test_spring_compress_paired_end {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true), 
        file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
    ]

    SPRING_COMPRESS ( input )
}
