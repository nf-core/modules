#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR       } from '../../../../../modules/nf-core/untar/main.nf'
include { KMCP_COMPUTE } from '../../../../../modules/nf-core/kmcp/compute/main.nf'
include { KMCP_INDEX } from '../../../../../modules/nf-core/kmcp/index/main.nf'
include { KMCP_SEARCH } from '../../../../../modules/nf-core/kmcp/search/main.nf'
include { KMCP_PROFILE } from '../../../../../modules/nf-core/kmcp/profile/main.nf'

workflow test_kmcp_profile {

    input_compute= UNTAR ( [
            [ id:'test' ],
            file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/kmcp/kmcp_profile.tar.gz", checkIfExists: true)
        ]).untar


    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    mode = 3

    KMCP_COMPUTE (input_compute )
    KMCP_INDEX ( KMCP_COMPUTE.out.outdir )
    KMCP_SEARCH ( KMCP_INDEX.out.kmcp.map{it[1]}, input )
    KMCP_PROFILE ( KMCP_SEARCH.out.result, input_compute.map{ it[1] }, mode )
}
