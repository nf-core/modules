#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GOLEFT_INDEXCOV } from '../../../../../modules/nf-core/goleft/indexcov/main.nf'

workflow test_goleft_indexcov {
    
   input = Channel.of([
        [ id:'test' ], // meta map
	    [
            file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true)
            ],
            [
            file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam_bai'], checkIfExists: true)
            ]
    ])

    fai = Channel.of([ [:], file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)])

    GOLEFT_INDEXCOV ( input , fai)
}
