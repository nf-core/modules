#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ILASTIK_MULTICUT } from '../../../../../modules/nf-core/ilastik/multicut/main.nf'

workflow test_ilastik_multicut {

    input = [
        [ id:'probabilities' ], // meta map
        file("/Users/florian_wuennemann/1_Projects/nf_core/test_data/wga.h5")
    ]
    ilp = [file("/Users/florian_wuennemann/1_Projects/nf_core/test_data/multicut_140rc8_plant.ilp")]
    probs = [file("/Users/florian_wuennemann/1_Projects/nf_core/test_data/wga-data_Probabilities.h5")]

    ILASTIK_MULTICUT ( input, ilp, probs )
}
