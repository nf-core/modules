#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR       } from '../../../../../modules/nf-core/untar/main.nf'
include { KMCP_COMPUTE } from '../../../../../modules/nf-core/kmcp/compute/main.nf'

workflow test_kmcp_compute {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    KMCP_COMPUTE ( input )
}

workflow test_kmcp_compute_directory {

input = UNTAR ( [
            [ id:'test' ],
            file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/kmcp/kmcp_profile.tar.gz", checkIfExists: true)
        ]).untar

KMCP_COMPUTE ( input )

}
