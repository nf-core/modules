#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ODGI_UNCHOP } from '../../../../../modules/nf-core/odgi/unchop/main.nf'

workflow test_odgi_unchop {

    input = [ [ id:'test' ],
              file(params.test_data['homo_sapiens']['pangenome']['pangenome_smoothxg_gfa'], checkIfExists: true),
            ]

    ODGI_UNCHOP ( input )
}
