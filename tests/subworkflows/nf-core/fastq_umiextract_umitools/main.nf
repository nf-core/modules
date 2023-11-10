#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQ_UMIEXTRACT_UMITOOLS } from '../../../../subworkflows/nf-core/fastq_umiextract_umitools/main.nf'

workflow test_fastq_umiextract_umitools {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)

        ]
    ]

    FASTQ_UMIEXTRACT_UMITOOLS ( input )
}

workflow test_fastq_umiextract_umitools_singleend {

    input = [
        [ id:'test', single_end:true ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)

        ]
    ]

    FASTQ_UMIEXTRACT_UMITOOLS ( input )
}
