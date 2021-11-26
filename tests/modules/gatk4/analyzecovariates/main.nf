#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_ANALYZECOVARIATES } from '../../../../modules/gatk4/analyzecovariates/main.nf'

workflow test_gatk4_analyzecovariates {
    input     = [ [ id:'test' ], // meta map
                    file(params.test_data['homo_sapiens']['illumina']['test_baserecalibrator_table'], checkIfExists: true),
                    [],
                    []
                  ]

    csvout = true
    ignoretimewarning = true

    GATK4_ANALYZECOVARIATES ( input, [], [] )
}
