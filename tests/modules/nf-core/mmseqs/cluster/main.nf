#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR          } from '../../../../../modules/nf-core/untar/main.nf'
include { MMSEQS_CLUSTER } from '../../../../../modules/nf-core/mmseqs/cluster/main.nf'

workflow test_mmseqs_cluster {

    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['genome']['mmseqs_tar_gz'], checkIfExists: true)
    ]

    ch_db = UNTAR ( input ).untar

    MMSEQS_CLUSTER ( ch_db )

}
