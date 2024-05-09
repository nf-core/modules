#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { GETORGANELLE_CONFIG }    from '../../../../../modules/nf-core/getorganelle/config/main.nf'
include { GETORGANELLE_FROMREADS } from '../../../../../modules/nf-core/getorganelle/fromreads/main.nf'

workflow test_getorganelle_fromreads {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    db = Channel.of('animal_mt')

    GETORGANELLE_CONFIG    ( db )
    GETORGANELLE_FROMREADS ( input , GETORGANELLE_CONFIG.out.db)
}
