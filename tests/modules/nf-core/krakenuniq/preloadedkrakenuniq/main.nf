#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KRAKENUNIQ_PRELOADEDKRAKENUNIQ } from '../../../../../modules/nf-core/krakenuniq/preloadedkrakenuniq/main.nf'

workflow test_krakenuniq_preloadedkrakenuniq {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    KRAKENUNIQ_PRELOADEDKRAKENUNIQ ( input )
}
