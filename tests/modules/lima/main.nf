#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { LIMA } from '../../../modules/lima/main.nf' addParams( options: [args: '--isoseq --peek-guess', suffix: ".fl"] )

workflow test_lima {

    input = [
                [ id:'test' ], // meta map
                file(params.test_data['homo_sapiens']['pacbio']['ccs'],     checkIfExists: true),
            ]
    primers = [ file(params.test_data['homo_sapiens']['pacbio']['primers'], checkIfExists: true) ]

    LIMA ( input, primers )
}
