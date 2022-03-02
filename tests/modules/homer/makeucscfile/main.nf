#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HOMER_MAKETAGDIRECTORY } from '../../../../modules/homer/maketagdirectory/main.nf'
include { HOMER_MAKEUCSCFILE } from '../../../../modules/homer/makeucscfile/main.nf'

workflow test_homer_makeucscfile {
    input = [[id:'test'],
             [file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true),
              file(params.test_data['sarscov2']['genome']['test2_bed'], checkIfExists: true)]]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    HOMER_MAKETAGDIRECTORY (input, fasta)
    HOMER_MAKEUCSCFILE ( HOMER_MAKETAGDIRECTORY.out.tagdir )
}

