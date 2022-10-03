#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTK_FASTK } from '../../../../../modules/nf-core/fastk/fastk/main.nf'
include { FASTK_MERGE } from '../../../../../modules/nf-core/fastk/merge/main.nf'

workflow test_fastk_merge_hist_only {

    input1 = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    input2= [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    FASTK_FASTK ( Channel.of( input1, input2 ) )
    FASTK_MERGE (
        FASTK_FASTK.out.hist.groupTuple()
            .join( FASTK_FASTK.out.ktab.groupTuple(), remainder: true )
            .join( FASTK_FASTK.out.prof.groupTuple(), remainder: true )
            .map { meta, hist, ktab, prof -> [meta, hist, ktab ? ktab.flatten() : [] , prof ? prof.flatten() : [] ] }
    )
}

workflow test_fastk_merge_all_files {

    input1 = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    input2= [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    FASTK_FASTK ( Channel.of( input1, input2 ) )
    FASTK_MERGE (
        FASTK_FASTK.out.hist.groupTuple()
            .join( FASTK_FASTK.out.ktab.groupTuple(), remainder: true )
            .join( FASTK_FASTK.out.prof.groupTuple(), remainder: true )
            .map { meta, hist, ktab, prof -> [meta, hist, ktab ? ktab.flatten() : [] , prof ? prof.flatten() : [] ] }
    )
}
