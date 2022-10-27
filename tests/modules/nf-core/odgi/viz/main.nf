#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { ODGI_VIZ } from "$moduleDir/modules/nf-core/odgi/viz/main.nf"

workflow test_odgi_viz {

    input = [ [ id:'test' ],
              [ file(params.test_data['homo_sapiens']['pangenome']['pangenome_smoothxg_gfa'], checkIfExists: true) ],
            ]

    ODGI_VIZ ( input )
}
