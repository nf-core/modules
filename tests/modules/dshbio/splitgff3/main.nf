#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DSHBIO_SPLITGFF3 } from '../../../../modules/dshbio/splitgff3/main.nf'

workflow test_dshbio_splitgff3 {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['genome']['genome_gff3_gz'], checkIfExists: true) ]
            ]

    DSHBIO_SPLITGFF3 ( input )
}
