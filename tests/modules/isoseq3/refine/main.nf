#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ISOSEQ3_REFINE } from '../../../../modules/isoseq3/refine/main'

workflow test_isoseq3_refine {

    input = [
            [ id:'test' ], // meta map
            file(params.test_data['homo_sapiens']['pacbio']['lima'], checkIfExists: true),
    ]
    primers = file(params.test_data['homo_sapiens']['pacbio']['primers'], checkIfExists: true)

    ISOSEQ3_REFINE ( input, primers )
}
