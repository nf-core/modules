#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ODGI_DRAW } from '../../../../../modules/nf-core/odgi/draw/main.nf'

workflow test_odgi_draw {

    input = [ [ id:'test' ],
              file(params.test_data['homo_sapiens']['pangenome']['odgi']['pangenome_og'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['pangenome']['odgi']['pangenome_lay'], checkIfExists: true)
            ]

    ODGI_DRAW ( input )
}
