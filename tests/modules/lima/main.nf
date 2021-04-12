#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { LIMA } from '../../../software/lima/main.nf' addParams( options: [:] )

workflow test_lima {

    input = [ [ id:'test' ], // meta map
              file(params.test_data['homo_sapiens']['pacbio']['ccs'], checkIfExists: true)
            ]
    primers = file(params.test_data['homo_sapiens']['pacbio']['primers'], checkIfExists: true)

    LIMA ( input, primers )
}
