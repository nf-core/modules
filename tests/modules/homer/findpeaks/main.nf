#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HOMER_MAKETAGDIRECTORY } from '../../../../modules/homer/maketagdirectory/main.nf' addParams( options: [args: '-format bed'] )
include { HOMER_FINDPEAKS } from '../../../../modules/homer/findpeaks/main.nf' addParams( options: [args: '-style factor'] )

workflow test_homer_findpeaks {
    input = [[id:'test'],
             [file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true),
              file(params.test_data['sarscov2']['genome']['test2_bed'], checkIfExists: true)]]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    HOMER_MAKETAGDIRECTORY (input, fasta)
    HOMER_FINDPEAKS ( HOMER_MAKETAGDIRECTORY.out.tagdir )
}

