#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MCMICRO_SCIMAP } from '../../../../../modules/nf-core/mcmicro/scimap/main.nf'

workflow test_mcmicro_scimap {

    input = [
        [ id:'test', single_end:false ], // meta map
        file("/Users/luiskuhn/Downloads/cycif_tonsil_small--unmicst_cell.csv", checkIfExists: true)
    ]

    MCMICRO_SCIMAP ( input )
}
