#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MSISENSOR2_MSI } from '../../../../modules/msisensor2/msi/main.nf'

workflow test_msisensor2_msi {

    input = [// meta map
        [ id:'test'],
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_cram'], checkIfExists: true),
        []
    ]

    MSISENSOR2_MSI ( input )
}
