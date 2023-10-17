#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FAMSA_ALIGN     } from '../../../../../modules/nf-core/famsa/align/main.nf'
include { TCOFFEE_ALNCOMPARE } from '../../../../../modules/nf-core/tcoffee/alncompare/main.nf'

workflow test_tcoffee_alncompare {

    input = [
        [ id:'test' ], // meta map
        file("https://raw.githubusercontent.com/nf-core/test-datasets/multiplesequencealign/testdata/setoxin-ref.fa", checkIfExists: true)
    ]
    input2= [
        [ id:'test' ], // meta map
        file("https://github.com/nf-core/test-datasets/blob/multiplesequencealign/testdata/setoxin.ref", checkIfExists: true)
    ]

    ch_aln=FAMSA_ALIGN ( input , [[:],[]] ).alignment
    ch_to_eval=ch_aln.combine(input2, by:0)


}

