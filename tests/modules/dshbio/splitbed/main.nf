#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DSHBIO_SPLITBED } from '../../../../modules/dshbio/splitbed/main.nf'

workflow test_dshbio_splitbed {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true) ]
            ]

    DSHBIO_SPLITBED ( input )
}
