#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HOMER_MAKETAGDIRECTORY } from '../../../../software/homer/maketagdirectory/main.nf' addParams( options: [:] )
include { HOMER_MAKEUCSCFILE } from '../../../../software/homer/makeucscfile/main.nf' addParams( options: [:] )

workflow test_homer_makeucscfile {
    def input = []
    input = [[id: 'test'], // meta map
             [file("${launchDir}/tests/data/bed/A.bed", checkIfExists: true),
              file("${launchDir}/tests/data/bed/B.bed", checkIfExists: true)]]

    fasta = file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true)

    HOMER_MAKETAGDIRECTORY (input, fasta)
    HOMER_MAKEUCSCFILE ( HOMER_MAKETAGDIRECTORY.out.tagdir )
}

