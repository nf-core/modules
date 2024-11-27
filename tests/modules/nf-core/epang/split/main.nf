#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MAFFT_ALIGN as MAFFT } from '../../../../../modules/nf-core/mafft/align/main.nf'
include { EPANG_SPLIT          } from '../../../../../modules/nf-core/epang/split/main.nf'

workflow test_epang_split {

    input = [
        [ id:'test', single_end:false ], // meta map
        file('https://raw.githubusercontent.com/nf-core/test-datasets/phyloplace/testdata/PF14720_seed.alnfaa', checkIfExists: true)
    ]

    MAFFT (
        input,
        [[id:"test2"], file('https://raw.githubusercontent.com/nf-core/test-datasets/phyloplace/testdata/PF14720_3_sequences.faa', checkIfExists: true)],
        [[:],[]],
        [[:],[]],
        [[:],[]],
        [[:],[]],
        false
    )

    EPANG_SPLIT (
        MAFFT.out.fas.map {
            [
                [ id:'test'],
                file('https://raw.githubusercontent.com/nf-core/test-datasets/phyloplace/testdata/PF14720_seed.alnfaa', checkIfExists: true),
                it[1]
            ]
        }
    )
}
