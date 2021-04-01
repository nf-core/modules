#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DSH_FILTERBED } from '../../../../software/dsh/filterbed/main.nf' addParams( options: [suffix: '.filtered', args: '--range chr1:0-1000'] )

workflow test_dsh_filterbed {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true) ]
            ]

    DSH_FILTERBED ( input )
}
