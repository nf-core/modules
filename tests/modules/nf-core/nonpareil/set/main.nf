#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NONPAREIL_NONPAREIL } from '../../../../../modules/nf-core/nonpareil/nonpareil/main.nf'
include { NONPAREIL_SET } from '../../../../../modules/nf-core/nonpareil/set/main.nf'

workflow test_nonpareil_set {

    input = Channel.fromList([
            [[ id:'test', single_end:false ], file(params.test_data['candidatus_portiera_aleyrodidarum']['illumina']['test_1_fastq_gz'], checkIfExists: true)],
            [[ id:'test2', single_end:false ], file(params.test_data['candidatus_portiera_aleyrodidarum']['illumina']['test_1_fastq_gz'], checkIfExists: true)],
        ]
    )

    NONPAREIL_NONPAREIL ( input, 'fastq', 'kmer' )
    NONPAREIL_SET (NONPAREIL_NONPAREIL.out.npo.map{meta, npo -> [[id:'test'], npo]}.groupTuple())

}
