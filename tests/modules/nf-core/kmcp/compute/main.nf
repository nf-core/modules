#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KMCP_COMPUTE } from '../../../../../modules/nf-core/kmcp/compute/main.nf'

workflow test_kmcp_compute {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    KMCP_COMPUTE ( input, [], [] )
}
