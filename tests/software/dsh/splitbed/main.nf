#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DSH_SPLITBED } from '../../../../software/dsh/splitbed/main.nf' addParams( options: [suffix: '.', args: '--records 2'] )

workflow test_dsh_splitbed {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true) ]
            ]

    DSH_SPLITBED ( input )
}
