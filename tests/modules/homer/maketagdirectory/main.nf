#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HOMER_MAKETAGDIRECTORY } from '../../../../modules/homer/maketagdirectory/main.nf' addParams( options: [:] )

workflow test_homer_maketagdirectory {
    input = [[id:'test'],
             [file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)]]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    HOMER_MAKETAGDIRECTORY (input, fasta)
}


workflow test_homer_meta_maketagdirectory {
    input =
        [[[ id:'test1'],
          [file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true)]],
         [[ id:'test2'],
          [file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)]]]

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    meta_input = [[id: 'meta_test']] + [ input.collect{it[1]}.flatten() ]

    HOMER_MAKETAGDIRECTORY (meta_input, fasta)
}
