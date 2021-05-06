#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PBBAM_PBMERGE } from '../../../../software/pbbam/pbmerge/main.nf' addParams( options: [:] )

workflow test_pbbam_pbmerge {

    input = [
                [ id:'test' ], // meta map
                [
                    file(params.test_data['homo_sapiens']['pacbio']['alz1000'],   checkIfExists: true),
                    file(params.test_data['homo_sapiens']['pacbio']['alz100000'], checkIfExists: true)
                ]
            ]

    PBBAM_PBMERGE ( input )
}
