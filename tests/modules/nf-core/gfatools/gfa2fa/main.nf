#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GFATOOLS_GFA2FA } from '../../../../../modules/nf-core/gfatools/gfa2fa/main.nf'

workflow test_gfatools_gfa2fa {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['assembly_gfa'], checkIfExists: true)
    ]

    GFATOOLS_GFA2FA ( input )
}
