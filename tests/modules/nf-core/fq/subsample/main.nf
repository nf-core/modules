#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FQ_SUBSAMPLE                              } from '../../../../../modules/nf-core/fq/subsample/main.nf'
include { FQ_SUBSAMPLE as FQ_SUBSAMPLE_PROBABILITY  } from '../../../../../modules/nf-core/fq/subsample/main.nf'
include { FQ_SUBSAMPLE as FQ_SUBSAMPLE_RECORD_COUNT } from '../../../../../modules/nf-core/fq/subsample/main.nf'
include { FQ_SUBSAMPLE as FQ_SUBSAMPLE_NO_ARGS      } from '../../../../../modules/nf-core/fq/subsample/main.nf'

workflow test_fq_subsample_no_args {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    FQ_SUBSAMPLE_NO_ARGS( input )
}

workflow test_fq_subsample_probability {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    FQ_SUBSAMPLE_PROBABILITY ( input )
}

workflow test_fq_subsample_record_count {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    FQ_SUBSAMPLE_RECORD_COUNT ( input )
}

workflow test_fq_subsample_single {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
        ]
    ]
    FQ_SUBSAMPLE ( input )
}

workflow test_fq_subsample_no_input {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        []
    ]
    FQ_SUBSAMPLE ( input )
}