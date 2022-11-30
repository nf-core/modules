#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SORTMERNA } from '../../../../modules/nf-core/sortmerna/main.nf'

//
// Test with single-end data
//
workflow test_sortmerna_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]
    fasta = [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]
    SORTMERNA ( input, fasta )
}

//
// Test with paired-end data
//
workflow test_sortmerna_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
    fasta = [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]
    SORTMERNA ( input, fasta )
}

