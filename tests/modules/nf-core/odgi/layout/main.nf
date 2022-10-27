#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { ODGI_LAYOUT } from "$moduleDir/modules/nf-core/odgi/layout/main.nf"

workflow test_odgi_layout {

    input = [ [ id:'test' ],
              file(params.test_data['homo_sapiens']['pangenome']['pangenome_smoothxg_gfa'], checkIfExists: true),
            ]

    ODGI_LAYOUT ( input )
}
