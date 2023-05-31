#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FALCO } from '../../../../modules/nf-core/falco/main.nf'

//
// Test with single-end data
//
workflow test_falco_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]

    FALCO ( input )
}

//
// Test with paired-end data
//
workflow test_falco_paired_end {
    input = [ [id: 'test', single_end: false], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]

    FALCO ( input )
}

workflow test_falco_interleaved {
    input = [
                [id: 'test', single_end: false], // meta map
                [file(params.test_data['sarscov2']['illumina']['test_interleaved_fastq_gz'], checkIfExists: true)]
            ]

    FALCO ( input )
}
