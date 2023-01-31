#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VG_INDEX } from '../../../../../modules/nf-core/vg/index/main.nf'

workflow test_vg_index {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['pangenome']['pangenome_gfaffix_gfa'], checkIfExists: true),
        []
    ]

    mapping = [[],[]]

    VG_INDEX (
        input,
        mapping
    )
}
