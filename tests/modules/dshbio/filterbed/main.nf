#!/usr/bin/env nextflow



include { DSHBIO_FILTERBED } from '../../../../modules/dshbio/filterbed/main.nf'

workflow test_dshbio_filterbed {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true) ]
            ]

    DSHBIO_FILTERBED ( input )
}
