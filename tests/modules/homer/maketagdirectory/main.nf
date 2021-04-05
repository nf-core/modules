#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HOMER_MAKETAGDIRECTORY } from '../../../../software/homer/maketagdirectory/main.nf' addParams( options: [:] )

workflow test_homer_maketagdirectory {

    def input = []
    input = [[id: 'test'], // meta map
             [file("${launchDir}/tests/data/bed/A.bed", checkIfExists: true),
              file("${launchDir}/tests/data/bed/B.bed", checkIfExists: true)]]

    fasta = file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true)

    HOMER_MAKETAGDIRECTORY (input, fasta)
}


workflow test_homer_meta_maketagdirectory {

    def input = []
    input = [[[id: 'test1'], // meta map
              [file("${launchDir}/tests/data/bed/A.bed", checkIfExists: true)]],
             [[id: 'test2'],
              [file("${launchDir}/tests/data/bed/B.bed", checkIfExists: true)]]]

    fasta = file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true)

    meta_input = [[id: 'meta_test']] + [ input.collect{it[1]}.flatten() ]

    HOMER_MAKETAGDIRECTORY (meta_input, fasta)
}
