#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HICEXPLORER_HICPCA } from '../../../../../modules/nf-core/hicexplorer/hicpca/main.nf'

workflow test_hicexplorer_hicpca {

    input = [
        [ id:'test' ], // meta map
        file('https://github.com/deeptools/HiCExplorer/raw/master/hicexplorer/test/test_data/hicPCA/mm9_reduced_chr1.cool', checkIfExists: true)
    ]

    HICEXPLORER_HICPCA ( input )
}
