#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SRAHUMANSCRUBBER_INITDB } from '../../../../../modules/nf-core/srahumanscrubber/initdb/main.nf'
include { SRAHUMANSCRUBBER_SCRUB } from '../../../../../modules/nf-core/srahumanscrubber/scrub/main.nf'

workflow test_srahumanscrubber_scrub {
    
    input = [
        [id: 'test', single_end: false], // meta map
        [file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
         file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)]
    ]

    SRAHUMANSCRUBBER_INITDB()
    SRAHUMANSCRUBBER_SCRUB ( input, SRAHUMANSCRUBBER_INITDB.out.db )
}
