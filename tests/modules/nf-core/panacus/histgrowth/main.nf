#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PANACUS_HISTGROWTH } from '../../../../../modules/nf-core/panacus/histgrowth/main.nf'

workflow test_panacus_histgrowth {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['pangenome']['pangenome_smoothxg_gfa'], checkIfExists: true),
    ]

    PANACUS_HISTGROWTH ( input, [], [], [] )
}
