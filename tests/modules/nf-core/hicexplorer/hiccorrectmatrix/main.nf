#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HICEXPLORER_HICCORRECTMATRIX } from '../../../../../modules/nf-core/hicexplorer/hiccorrectmatrix/main.nf'
include { COOLER_DUMP    } from '../../../../../modules/nf-core/cooler/dump/main.nf'

workflow test_hicexplorer_hiccorrectmatrix {

    input = [
        [ id:'test', single_end:false ], // meta map
        file('https://github.com/deeptools/HiCExplorer/raw/master/hicexplorer/test/test_data/hicCorrectMatrix/gm12878_raw_values.cool', checkIfExists: true)
    ]

    HICEXPLORER_HICCORRECTMATRIX ( input )
    COOLER_DUMP ( HICEXPLORER_HICCORRECTMATRIX.out.corrected.combine([[:]]) )
}
