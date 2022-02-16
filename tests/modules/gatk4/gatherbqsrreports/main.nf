#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_GATHERBQSRREPORTS } from '../../../../modules/gatk4/gatherbqsrreports/main.nf'

workflow test_gatk4_gatherbqsrreports {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_baserecalibrator_table'], checkIfExists: true)
    ]

    GATK4_GATHERBQSRREPORTS ( input )
}

workflow test_gatk4_gatherbqsrreports_multiple {

    input = [
        [ id:'test', single_end:false ], // meta map
        [file(params.test_data['homo_sapiens']['illumina']['test_baserecalibrator_table'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_baserecalibrator_table'], checkIfExists: true)
        ]
    ]

    GATK4_GATHERBQSRREPORTS ( input )
}
