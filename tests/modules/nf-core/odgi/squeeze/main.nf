#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ODGI_SQUEEZE } from '../../../../../modules/nf-core/odgi/squeeze/main.nf'

workflow test_odgi_squeeze {
    odgi = file(params.test_data['homo_sapiens']['pangenome']['odgi']['pangenome_og'], checkIfExists: true)

    input = [ [ id:'test' ], // meta map
              [ odgi ]
            ]

    ODGI_SQUEEZE ( input )
}
