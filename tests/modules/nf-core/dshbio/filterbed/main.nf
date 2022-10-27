#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { DSHBIO_FILTERBED } from "$moduleDir/modules/nf-core/dshbio/filterbed/main.nf"

workflow test_dshbio_filterbed {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true) ]
            ]

    DSHBIO_FILTERBED ( input )
}
