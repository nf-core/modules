#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ODGI_SQUEEZE } from '../../../../../modules/nf-core/odgi/squeeze/main.nf'

workflow test_odgi_squeeze {
    odgi = file(params.test_data['homo_sapiens']['pangenome']['odgi']['pangenome_og'], checkIfExists: true)
    gfa = file(params.test_data['homo_sapiens']['pangenome']['pangenome_smoothxg_gfa'], checkIfExists: true)

    input = [ [ id:'test' ], // meta map
              [ odgi,
                gfa]
            ]

    ODGI_SQUEEZE ( input )
}
