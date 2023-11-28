#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQTK_TRIM as SEQTK_TRIM_SINGLE } from '../../../../../modules/nf-core/seqtk/trim/main.nf'
include { SEQTK_TRIM as SEQTK_TRIM_PAIRED } from '../../../../../modules/nf-core/seqtk/trim/main.nf'


//
// Test with single-end data
//
workflow test_seqtk_trim_single_end {

    input = [ 
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    SEQTK_TRIM_SINGLE ( input )
}

//
// Test with paired-end data
//
workflow test_seqtk_trim_paired_end {

    input = [ 
        [ id:'test', single_end:false ], // meta map
        [ 
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) 
        ]
    ]

    SEQTK_TRIM_PAIRED ( input )
}