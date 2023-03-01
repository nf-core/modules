#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HICEXPLORER_HICFINDTADS } from '../../../../../modules/nf-core/hicexplorer/hicfindtads/main.nf'

workflow test_hicexplorer_hicfindtads {

    input = [
        [ id:'test', single_end:false ], // meta map
        file('https://github.com/deeptools/HiCExplorer/raw/master/hicexplorer/test/test_data/hicTADClassifier/gm12878_chr1.cool', checkIfExists: true)
    ]

    HICEXPLORER_HICFINDTADS ( input , 20000)
}
