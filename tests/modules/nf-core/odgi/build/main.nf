#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ODGI_BUILD } from '../../../../../modules/nf-core/odgi/build/main.nf'

workflow test_odgi_build {

    input = [ [ id:'test' ],
              file(params.test_data['homo_sapiens']['pangenome']['pangenome_smoothxg_gfa'], checkIfExists: true),
            ]

    ODGI_BUILD ( input )
}
