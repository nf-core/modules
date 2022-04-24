#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MOTUS_PROFILE } from '../../../../modules/motus/profile/main.nf'

workflow test_motus_profile_single_end {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    MOTUS_PROFILE ( input )
}

workflow test_motus_profile_paired_end {

    input = [
        [ id:'test', single_end:false ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
          file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
    ]
    MOTUS_PROFILE ( input )

}
