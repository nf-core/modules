#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQTK_MERGEPE                                           } from '../../../../modules/seqtk/mergepe/main.nf'
include { KHMER_NORMALIZEBYMEDIAN                                 } from '../../../../modules/khmer/normalizebymedian/main.nf'
include { KHMER_NORMALIZEBYMEDIAN as KHMER_NORMALIZEBYMEDIAN_ARGS } from '../../../../modules/khmer/normalizebymedian/main.nf'

workflow test_khmer_normalizebymedian_only_pe {
    
    pe_reads = [
        [ id:'khmer_test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) 
        ]
    ]

    SEQTK_MERGEPE(pe_reads)

    KHMER_NORMALIZEBYMEDIAN ( SEQTK_MERGEPE.out.reads.collect { it[1] }, [], 'only_pe' )
}

workflow test_khmer_normalizebymedian_only_se {
    
    se_reads = [
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) 
    ]

    KHMER_NORMALIZEBYMEDIAN ( [], se_reads, 'only_se' )
}

workflow test_khmer_normalizebymedian_mixed {
    
    pe_reads = [
        [ id:'khmer_test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) 
        ]
    ]
    se_reads = file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)

    SEQTK_MERGEPE(pe_reads)

    KHMER_NORMALIZEBYMEDIAN ( SEQTK_MERGEPE.out.reads.map { it[1] }, se_reads, 'mixed' )
}

workflow test_khmer_normalizebymedian_multiple_pe {
    
    pe_reads = [
        [ id:'khmer_test0', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) 
        ],
        [ id:'khmer_test1', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) 
        ]
    ]
    se_reads = file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)

    SEQTK_MERGEPE(pe_reads)

    KHMER_NORMALIZEBYMEDIAN ( SEQTK_MERGEPE.out.reads.collect { it[1] }, se_reads, 'multiple_pe' )
}

workflow test_khmer_normalizebymedian_args {
    
    pe_reads = [
        [ id:'khmer_test0', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) 
        ]
    ]
    se_reads = file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)

    SEQTK_MERGEPE(pe_reads)

    KHMER_NORMALIZEBYMEDIAN_ARGS ( SEQTK_MERGEPE.out.reads.collect { it[1] }, se_reads, 'args' )
}
