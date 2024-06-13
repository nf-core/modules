#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GFATOOLS_STAT } from '../../../../../modules/nf-core/gfatools/stat/main.nf'

workflow test_gfatools_stat {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['assembly_gfa'], checkIfExists: true)
    ]

    GFATOOLS_STAT ( input )
}
