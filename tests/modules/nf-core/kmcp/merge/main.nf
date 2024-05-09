#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KMCP_COMPUTE } from '../../../../../modules/nf-core/kmcp/compute/main.nf'
include { KMCP_INDEX } from '../../../../../modules/nf-core/kmcp/index/main.nf'
include { KMCP_SEARCH } from '../../../../../modules/nf-core/kmcp/search/main.nf'
include { KMCP_MERGE } from '../../../../../modules/nf-core/kmcp/merge/main.nf'

workflow test_kmcp_merge {

    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    KMCP_COMPUTE ( input )
    KMCP_INDEX ( KMCP_COMPUTE.out.outdir )
    KMCP_SEARCH ( KMCP_INDEX.out.kmcp.map{it[1]}, input )
    KMCP_MERGE(KMCP_SEARCH.out.result)
}

