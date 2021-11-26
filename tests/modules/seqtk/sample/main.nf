#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQTK_SAMPLE } from '../../../../modules/seqtk/sample/main.nf'

//
// Test with single-end data
//
workflow test_seqtk_sample_single_end {

    input = [ [ id:'test', single_end:true ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]

    SEQTK_SAMPLE ( input, 50 )
}

//
// Test with paired-end data
//
workflow test_seqtk_sample_paired_end {

    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]

    SEQTK_SAMPLE ( input, 50 )
}
