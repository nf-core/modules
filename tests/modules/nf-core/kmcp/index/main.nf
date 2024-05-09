#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR       } from '../../../../../modules/nf-core/untar/main.nf'
include { KMCP_COMPUTE } from '../../../../../modules/nf-core/kmcp/compute/main.nf'
include { KMCP_INDEX   } from '../../../../../modules/nf-core/kmcp/index/main.nf'

workflow test_kmcp_index {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    KMCP_COMPUTE ( input )
    KMCP_INDEX ( KMCP_COMPUTE.out.outdir )
}

