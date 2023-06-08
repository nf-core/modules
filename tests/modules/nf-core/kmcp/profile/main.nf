#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR       } from '../../../../../modules/nf-core/untar/main.nf'
include { KMCP_COMPUTE } from '../../../../../modules/nf-core/kmcp/compute/main.nf'
include { KMCP_INDEX } from '../../../../../modules/nf-core/kmcp/index/main.nf'
include { KMCP_SEARCH } from '../../../../../modules/nf-core/kmcp/search/main.nf'
include { KMCP_PROFILE } from '../../../../../modules/nf-core/kmcp/profile/main.nf'

workflow test_kmcp_profile {

    input_compute = [
        [ id:'test', single_end:false ], // meta map
        file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/kmcp/NC_045512.2.fasta', checkIfExists: true)
    ]

    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    taxdump =  [ [],file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/kmcp/kmcp_profile.tar.gz', checkIfExists: true)]
    ch_taxid = file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/kmcp/seqid2taxid.map', checkIfExists: true)
    mode = 3

    UNTAR (taxdump)
    KMCP_COMPUTE ( input_compute )
    KMCP_INDEX ( KMCP_COMPUTE.out.outdir )
    KMCP_SEARCH ( KMCP_INDEX.out.kmcp.map{it[1]}, input )
    KMCP_PROFILE ( KMCP_SEARCH.out.result, UNTAR.out.untar.map{ it[1] }, ch_taxid, mode )
}
