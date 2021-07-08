#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HOMER_MAKETAGDIRECTORY } from '../../../../modules/homer/maketagdirectory/main.nf' addParams( options: [:] )
include { HOMER_FINDPEAKS } from '../../../../modules/homer/findpeaks/main.nf' addParams( options: [args: '-style factor'] )

workflow test_homer_findpeaks {
    input = [ [ id:'test'],
             file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true),
             file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    HOMER_MAKETAGDIRECTORY (input, fasta)
    HOMER_FINDPEAKS ( HOMER_MAKETAGDIRECTORY.out.tagdir )
}

