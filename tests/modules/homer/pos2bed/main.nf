#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HOMER_MAKETAGDIRECTORY } from '../../../../modules/homer/maketagdirectory/main.nf'
include { HOMER_FINDPEAKS } from '../../../../modules/homer/findpeaks/main.nf'
include { HOMER_POS2BED } from '../../../../modules/homer/pos2bed/main.nf'

workflow test_homer_pos2bed {
    input = [[id:'test'],
             [file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true),
              file(params.test_data['sarscov2']['genome']['test2_bed'], checkIfExists: true)]]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    HOMER_MAKETAGDIRECTORY (input, fasta)
    HOMER_FINDPEAKS ( HOMER_MAKETAGDIRECTORY.out.tagdir )

    HOMER_POS2BED ( HOMER_FINDPEAKS.out.txt )
}
