#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DSHBIO_FILTERBED } from '../../../../modules/dshbio/filterbed/main.nf' addParams( options: [suffix: '.filtered', args: '--range chr1:0-1000'] )

workflow test_dshbio_filterbed {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true) ]
            ]

    DSHBIO_FILTERBED ( input )
}
