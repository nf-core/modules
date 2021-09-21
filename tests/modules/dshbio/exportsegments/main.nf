#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DSHBIO_EXPORTSEGMENTS } from '../../../../modules/dshbio/exportsegments/main.nf' addParams( options: [:] )

workflow test_dshbio_exportsegments {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['assembly_gfa'], checkIfExists: true) ]
            ]

    DSHBIO_EXPORTSEGMENTS ( input )
}
