#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_FASTQTOSAM } from '../../../../../modules/nf-core/picard/fastqtosam/main.nf'

workflow test_picard_fastqtosam_single {

    input = [
        [ id:'test', single_end:true ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
        ]
    ]

    PICARD_FASTQTOSAM ( input )
}

workflow test_picard_fastqtosam_paired {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    PICARD_FASTQTOSAM ( input )
}

workflow test_picard_fastqtosam_paired_custom_samplename {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    PICARD_FASTQTOSAM ( input )
}
