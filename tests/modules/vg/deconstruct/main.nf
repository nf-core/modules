#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VG_DECONSTRUCT } from '../../../../modules/vg/deconstruct/main.nf'

workflow test_vg_deconstruct {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['homo_sapiens']['pangenome']['pangenome_fa'], checkIfExists: true)],
              [ file(params.test_data['homo_sapiens']['pangenome']['pangenome_seqwish_gfa'], checkIfExists: true) ]
            ]

    VG_DECONSTRUCT( input )
}
