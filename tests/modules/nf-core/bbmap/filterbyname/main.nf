#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BBMAP_FILTERBYNAME as BBMAP_FILTERBYNAME } from '../../../../../modules/nf-core/bbmap/filterbyname/main.nf'

workflow test_bbmap_filterbyname_paired_end {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    BBMAP_FILTERBYNAME ( input )
}

workflow test_bbmap_filterbyname_paired_end_interleaved_filter_seqs {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    BBMAP_FILTERBYNAME ( input )
}

workflow test_bbmap_filterbyname_single_end {

    input = [
        [ id:'test', single_end:true ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
        ]
    ]

    BBMAP_FILTERBYNAME ( input )
}

workflow test_bbmap_filterbyname_single_end_filter_seqs {

    input = [
        [ id:'test', single_end:true ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
        ]
    ]

    BBMAP_FILTERBYNAME ( input )
}
