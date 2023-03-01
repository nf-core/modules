#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HICEXPLORER_HICPLOTMATRIX } from '../../../../../modules/nf-core/hicexplorer/hicplotmatrix/main.nf'

workflow test_hicexplorer_hicplotmatrix {

    input = [
        [ id:'test', single_end:false ], // meta map
        file('https://github.com/deeptools/HiCExplorer/raw/master/hicexplorer/test/test_data/Li_et_al_2015.h5', checkIfExists: true),
        []
    ]

    HICEXPLORER_HICPLOTMATRIX ( input )
}
