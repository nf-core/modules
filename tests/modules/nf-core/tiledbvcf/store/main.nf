#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TILEDBVCF_STORE } from '../../../../../modules/nf-core/tiledbvcf/store/main.nf'
include { TILEDBVCF_CREATE } from '../../../../../modules/nf-core/tiledbvcf/create/main.nf'

workflow test_tiledbvcf_store {

    vcf_tbi = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true)
    ]


    ch_input = [
        [id: 'test_datastore'],
        'my_vcf_dataset'
    ]

    TILEDBVCF_CREATE ( ch_input )

    ch_uri = TILEDBVCF_CREATE.out.uri

    TILEDBVCF_STORE ( vcf_tbi, ch_uri )
}
