#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { PBBAM_PBMERGE } from "$moduleDir/modules/nf-core/pbbam/pbmerge/main.nf"

workflow test_pbbam_pbmerge {

    input = [
        [ id:'test' ], // meta map
        [
            file(params.test_data['homo_sapiens']['pacbio']['cluster']   , checkIfExists: true),
            file(params.test_data['homo_sapiens']['pacbio']['singletons'], checkIfExists: true)
        ]
    ]

    PBBAM_PBMERGE ( input )
}
