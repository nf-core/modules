#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_REVERTSAM } from '../../../../../modules/nf-core/gatk4/revertsam/main.nf'

workflow test_gatk4_revertsam {
    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
            ]

    GATK4_REVERTSAM ( input )
}

workflow test_gatk4_revertsam_stubs {
    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
            ]

    GATK4_REVERTSAM ( input )
}
