#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TBPROFILER_PROFILE } from '../../../../modules/tbprofiler/profile/main.nf'

workflow test_tbprofiler_profile_illumina {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    TBPROFILER_PROFILE ( input )
}

workflow test_tbprofiler_profile_nanopore {

    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['nanopore']['test_fastq_gz'], checkIfExists: true)
    ]

    TBPROFILER_PROFILE ( input )
}
