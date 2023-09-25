#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GANON_BUILDCUSTOM } from '../../../../../modules/nf-core/ganon/buildcustom/main.nf'

workflow test_ganon_buildcustom {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
    ]

    GANON_BUILDCUSTOM ( input, [], [] )
}


