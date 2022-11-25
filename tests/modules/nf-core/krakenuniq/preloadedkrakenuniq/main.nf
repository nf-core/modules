#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KRAKENUNIQ_PRELOADEDKRAKENUNIQ } from '../../../../../modules/nf-core/krakenuniq/preloadedkrakenuniq/main.nf'

workflow test_krakenuniq_preloadedkrakenuniq_paired {

    input = [
        [id: 'test', single_end: false], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    db = []
    ram_chunk_size = '8GB'
    save_output_fastqs = true
    report_file = true
    save_output = true

    KRAKENUNIQ_PRELOADEDKRAKENUNIQ ( input, db, ram_chunk_size, save_output_fastqs, report_file, save_output )
}

workflow test_krakenuniq_preloadedkrakenuniq_single {

    input = [
        [id:'test', single_end:true], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    db = []
    ram_chunk_size = '8GB'
    save_output_fastqs = true
    report_file = true
    save_output = true

    KRAKENUNIQ_PRELOADEDKRAKENUNIQ ( input, db, ram_chunk_size, save_output_fastqs, report_file, save_output )
}
