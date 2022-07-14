#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    HOMER_MAKETAGDIRECTORY as HOMER_MAKETAGDIRECTORY_BED
    HOMER_MAKETAGDIRECTORY as HOMER_MAKETAGDIRECTORY_BAM
} from '../../../../modules/homer/maketagdirectory/main.nf'

workflow test_homer_maketagdirectory_bed {
    input = [[id:'test'],
             [file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true),
              file(params.test_data['sarscov2']['genome']['test2_bed'], checkIfExists: true)]]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    HOMER_MAKETAGDIRECTORY_BED (input, fasta)
}


workflow test_homer_maketagdirectory_meta {
    input =
        [[[ id:'test1'],
          [file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)]],
         [[ id:'test2'],
          [file(params.test_data['sarscov2']['genome']['test2_bed'], checkIfExists: true)]]]

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    meta_input = [[id: 'meta_test']] + [ input.collect{it[1]}.flatten() ]

    HOMER_MAKETAGDIRECTORY_BED (meta_input, fasta)
}

workflow test_homer_maketagdirectory_bam {
    input = [[id:'test'],
             [file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)]]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    HOMER_MAKETAGDIRECTORY_BAM (input, fasta)
}
