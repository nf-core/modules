#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DSHBIO_FILTERGFF3 } from '../../../../modules/dshbio/filtergff3/main.nf' addParams( options: [suffix: '.filtered', args: '--range MT192765.1:0-1000'] )

workflow test_dshbio_filtergff3 {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['genome']['genome_gff3_gz'], checkIfExists: true) ]
            ]

    DSHBIO_FILTERGFF3 ( input )
}
